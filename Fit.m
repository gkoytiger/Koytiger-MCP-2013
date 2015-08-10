function [info_Fits] = Fit(Cy5NormalizedData, AverageData_list, MaxConcentrations, PercentExclude) 
%Greg Koytiger June 2011

info_Fits(:,1:2) = unique(AverageData_list(:,1:2),'rows');
cell_Fits = cell(length(info_Fits),1);
cell_Gofs = cell(length(info_Fits),1);
SignalToNoise = zeros(length(info_Fits),1);

%ASSUMPTION Maximum possible kD is 1000 uM, further thresholding should be done;
Kd_threshold = 1000;

%Iterates through all the unique interacting pairs and fits a saturation binding curve
parfor interactionI = 1:length(info_Fits)
   
    concentration_Maximum = MaxConcentrations((MaxConcentrations(:,1) == info_Fits(interactionI,2)),2);  
    Xaxis = [.01 .1 .2 .5 1 2 3 5]' .* (concentration_Maximum / 5); %scales the Xaxis by real peptide concentration
    Yaxis  = Cy5NormalizedData(Cy5NormalizedData(:,1) == info_Fits(interactionI,1) & Cy5NormalizedData(:,2) == info_Fits(interactionI,2), 4);
    
    Xaxis = [0; Xaxis]; Yaxis=[Yaxis; 0;0;0;0]; %add 0,0 point needed to fit high affinity interactions
    
    %Quadruiplicates the Xaxis so it is the same size as the data axis
    Xaxis = sort([Xaxis; Xaxis; Xaxis; Xaxis], 'descend');

    
    Background = AverageData_list(AverageData_list(:,1) == info_Fits(interactionI,1) & AverageData_list(:,2) == info_Fits(interactionI,2),5);
       
    %Calculates a robust estimate of the fold over background of the titration
    SignalToNoise(interactionI) = trimmean(Background, PercentExclude);
    
    %Fits the titration to single site saturation binding equation
    [cell_Fits{interactionI} ,cell_Gofs{interactionI}] = fit_Kd(Xaxis, Yaxis, Kd_threshold);
    
end

%Second loop to build the data table (necessary to make first loop parallel)

for i = 1:length(info_Fits)
    beta = coeffvalues(cell_Fits{i});
    ci = confint(cell_Fits{i});
    info_Fits(i, 3) = beta(2);                            %Stores Kd value
    info_Fits(i, 4) = cell_Gofs{i}.rsquare;               %Stores R2 of fit
    info_Fits(i, 5) = SignalToNoise(i);                   %Stores Fold over Background of titration
    info_Fits(i, 6) = ci(1,2); info_Fits(i, 7) = ci(2,2); %Records the 95% CI of the Kd fit
end

end