function [RTK_p RTK_conn_onc RTK_conn_not_onc sorted_RTKs] = generateFigure2a(matrix_Kd, Peptide_list, RTK_list, Oncogene_list, Interaction_threshold)
%Greg Koytiger May 2012
%Finds the relationship between oncogenes and array connectivity
%Generates visualization of this relationship shown in figure 2a

%Isolates RTK peptides
RTKs = RTK_list(ismember(RTK_list(:,1), Peptide_list(:,1)),:);

%Splits affinities into 10 equivalent bins
bins = 0:(Interaction_threshold*1000)/10:Interaction_threshold*1000;

connectivity_Density = zeros(length(bins), length(RTKs(:,1)));

%Finds the number of binding partners for an RTK in each bin
for i = 2:length(bins)
    for j = 1:length(RTKs)
        RTK_Peptides = ismember(Peptide_list(:,1), RTKs(j,1));
        tempBinding = matrix_Kd(:,RTK_Peptides);
        connectivity_Density(i,j) = connectivity_Density(i-1, j) + sum(sum((tempBinding > bins(i-1) & tempBinding <= bins(i))));
    end
end

%Sorts RTKs by connectivity at 1uM
row = length(connectivity_Density(:,1));
[~, ix] = sort(connectivity_Density(row,:) ,'descend');

sorted_connectivity_Density = connectivity_Density(:,ix);
sorted_RTKs = RTKs(ix,:);
isOncogene = ismember(sorted_RTKs(:,1), Oncogene_list(:,1));

%Mann-Whitney U Test on whether the oncogene/non-oncogene groups are
%similarly connected
RTK_p = ranksum(sorted_connectivity_Density(row,isOncogene),sorted_connectivity_Density(row,~isOncogene), 0.05, 'method', 'exact');
RTK_conn_onc = median(sorted_connectivity_Density(row,isOncogene));
RTK_conn_not_onc = median(sorted_connectivity_Density(row,~isOncogene));

fig1 = plot(bins, sorted_connectivity_Density(:,1:40));

title('RTK Connectivity vs Interaction Threshold');
xlabel('Kd threshold (nM)');
ylabel('Number of Interactions');


legend(sorted_RTKs(:,1), 'Location','NorthEastOutside');
set(fig1(isOncogene), 'Color', [1 0 0]);
set(fig1(~isOncogene), 'Color', [.5 .5 .5]);
set(fig1,'LineWidth', 2)
set(gcf, 'Position', [0 50 1600 925], 'Resize', 'off','PaperPositionMode','auto');
print('-depsc2', 'figure/2/RTK_Conn_Onc');

%Makes boxplot to visualize the distributions
graph_label = cell(size(RTKs(:,1),1),1);
for i = 1:length(graph_label)
    if isOncogene(i)
        graph_label{i} = 'Oncogene';
    else
        graph_label{i} = 'Not Oncogene';
    end
end

figure; boxplot(sorted_connectivity_Density(row,:)',graph_label); title('RTK Connectivity'); xlabel({[num2str(RTK_conn_onc), ' vs ', num2str(RTK_conn_not_onc)] ; ['p value : ', num2str(RTK_p, '%4.2e')]});
print('-depsc2', 'figure/2/RTK_Conn_boxplot');

end
