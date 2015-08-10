function generate_figure_1a(Domain_list, PercentExclude, interaction_Kd, Data)
%Greg Koytiger June 2012
%Generates Figure 1 Panel A which shows an example of a peptide titration
%and fit for KIT-pY900
%Utilizes the same code for fitting as the script

PlateRow = 169;
Peptide = 125;


MaxConcentrations = Data{PlateRow,2};
ArrayType = Data{PlateRow,4};

RawData_list = Reformat(Data, PlateRow, ArrayType, PercentExclude);  
[AverageData_list Cy5NormalizedData] = DataAverage(RawData_list, PercentExclude); 

info_Fits(:,1:2) = unique(AverageData_list(:,1:2),'rows');
filter = ismember(info_Fits(:,1), cell2mat(Domain_list(:,1))) & ismember(info_Fits(:,2), Peptide);
info_Fits(~filter,:) = [];

cell_Fits = cell(length(info_Fits),1);
cell_Gofs = cell(length(info_Fits),1);
SignalToNoise = zeros(length(info_Fits),1);
fit_points = cell(length(info_Fits),1);


%ASSUMPTION Maximum possible kD is 1000 uM, further thresholding should be done;
Kd_threshold = 1000;

%Iterates through all the unique interacting pairs and fits a saturation binding curve

parfor interactionI = 1:length(info_Fits)
   
    concentration_Maximum = MaxConcentrations((MaxConcentrations(:,1) == info_Fits(interactionI,2)),2);  
    Xaxis = [.01 .1 .2 .5 1 2 3 5]' .* (concentration_Maximum / 5); %scales the Xaxis by real peptide concentration
    Yaxis  = Cy5NormalizedData(Cy5NormalizedData(:,1) == info_Fits(interactionI,1) & Cy5NormalizedData(:,2) == info_Fits(interactionI,2), 4);
    
    Xaxis = [0; Xaxis]; Yaxis=[Yaxis; 0;0;0;0]; %add 0,0 data point
    
    %Quadruiplicates the Xaxis so it is the same size as the data axis
    Xaxis = sort([Xaxis; Xaxis; Xaxis; Xaxis], 'descend');

    
    Background = AverageData_list(AverageData_list(:,1) == info_Fits(interactionI,1) & AverageData_list(:,2) == info_Fits(interactionI,2),5);
       
    %Calculates a robust estimate of the fold over background of the titration
    SignalToNoise(interactionI) = trimmean(Background, PercentExclude);
    
    %Fits the titration to single site saturation binding equation
    [cell_Fits{interactionI} ,cell_Gofs{interactionI}] = fit_Kd(Xaxis, Yaxis, Kd_threshold);
    fit_points{interactionI} = [Xaxis, Yaxis];
end

for i = 1:length(info_Fits)
    beta = coeffvalues(cell_Fits{i});
    ci = confint(cell_Fits{i});
    info_Fits(i, 3) = beta(2);                            %Stores Kd value
    info_Fits(i, 4) = cell_Gofs{i}.rsquare;               %Stores R2 of fit
    info_Fits(i, 5) = SignalToNoise(i);                   %Stores Fold over Background of titration
    info_Fits(i, 6) = ci(1,2); info_Fits(i, 7) = ci(2,2); %Records the 95% CI of the Kd fit
end


filter = SignalToNoise > 5  & info_Fits(:,4) > 0.90 & info_Fits(:,3) < interaction_Kd;

info_Fits(~filter,:) = []; fit_points(~filter) = []; cell_Fits(~filter) = [];

[~,idx]= sort(info_Fits(:,3));

info_Fits = info_Fits(idx,:); fit_points = fit_points(idx); cell_Fits = cell_Fits(idx);
figure;
lengthplot = ceil(size(info_Fits,1) / 7);

for i = 1:size(info_Fits,1)
        subplot(lengthplot, 7, i);

        X = unique(fit_points{i}(:,1));
        X = sort(X, 'descend');
        Y = X;
        error = X;
        k = 4;
        for j = 1:9      
            Y(j) = trimmean(fit_points{i}(k-3:k, 2), 25);
            error(j) = std(fit_points{i}(k-3:k, 2));
            k = k+4;
        end
        
        axis([0 5 0 max(Y) * 1.25]);
        fig1 = plot(cell_Fits{i}, X, Y);
        set(fig1, 'LineWidth', 3, 'Color', colorassign(info_Fits(i,3) * 1000));
       
        
        set(gca, 'YTick', []);
        
        legend('hide');
        ylabel('');
        
        hold on;
        errorbar(X, Y , error, '.');
        
        xlabel(['Kd = ', num2str(info_Fits(i,3) * 1000, 3), ' nM  R^2 = ', num2str(info_Fits(i,4),3)]);
        title(Domain_list(cell2mat(Domain_list(:,1)) == info_Fits(i,1), 3));
        axis([0 5 0 max(Y) * 1.25]);
        
        hold off;        
end
print('-depsc2', 'figure/1/KIT-pY900');

end