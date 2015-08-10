ANKS1A_454_info = info_Fits(ismember(info_Fits(:,1), cell2mat(Domain_list(:,1))) & info_Fits(:,2) == 2340,:);
ANKS1A_454_cell = cell_Fits(ismember(info_Fits(:,1), cell2mat(Domain_list(:,1))) & info_Fits(:,2) == 2340,:);
ANKS1A_454_gof = cell_Gofs(ismember(info_Fits(:,1), cell2mat(Domain_list(:,1))) & info_Fits(:,2) == 2340,:);
ANKS1A_454_SNR = SignalToNoise(ismember(info_Fits(:,1), cell2mat(Domain_list(:,1))) & info_Fits(:,2) == 2340,:);
ANKS1A_454_fitpoints = cell_Fitted_Points(ismember(info_Fits(:,1), cell2mat(Domain_list(:,1))) & info_Fits(:,2) == 2340,:);

close all;

for i = 1:length(ANKS1A_454_info)
    if(ANKS1A_454_gof{i}.rsquare > 0.93 && ANKS1A_454_SNR(i) > 5)
        figure;
        %plot(ANKS1A_454_fitpoints{i}{1}, ANKS1A_454_fitpoints{i}{2}, '.');
        X = unique(ANKS1A_454_fitpoints{i}{1});
        X = sort(X, 'descend');
        Y = X;
        error = X;
        k = 4;
        for j = 1:9      
            Y(j) = trimmean(ANKS1A_454_fitpoints{i}{2}(k-3:k), 25);
            error(j) = std(ANKS1A_454_fitpoints{i}{2}(k-3:k));
            k = k+4;
        end
        
        fig1 = plot(ANKS1A_454_cell{i});
        set(gca, 'YTick', []);
        set(fig1, 'LineWidth', 2);
        legend('hide');
        ylabel('');
        xlabel('');
        
        hold on;
        errorbar(X, Y , error, '.');
        
        title_label = strcat(Domain_list(cell2mat(Domain_list(:,1)) == ANKS1A_454_info(i,1), 2), '  Kd = ', num2str(ANKS1A_454_cell{i}.Kd * 1000, 3), ' nM    R^2 = ', num2str(ANKS1A_454_gof{i}.rsquare,2));
        title(title_label);
        axis([0 5 0 max(Y) * 1.25]);
        hold off;        
    end
end
    