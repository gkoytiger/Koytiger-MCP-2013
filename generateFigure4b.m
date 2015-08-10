function generateFigure4b(matrix_Kd, Peptide_list, Domain_list, Oncogene_list, Interaction_threshold)
%Makes figure 4b and analyzes the relationship between adaptor connectivity
%and oncogencity
%Note that the way statistics are done and rank is visualized is subtely
%different on account of ties, most notably with 0's. If we visualized tied
%ranks then everything will be bunched together but stastics are all
%performed accurately


%Eliminates tandem domains from analysis (to remove double counting effects
%of interactions)
isTandem = ~cellfun('isempty', strfind(Domain_list(:,3), '-NC'));
matrix_Kd(isTandem,:) = [];
Domain_list(isTandem,:) = [];
Interaction_threshold = Interaction_threshold * 1000; %Converts to nM

%Initializes variables
in_links = zeros(length(Domain_list),1);
out_links = zeros(length(Domain_list),1);
peptide_count = zeros(length(Domain_list),1);
normalized_out_links = zeros(length(Domain_list),1);

%Finds the number of in-links and out-links
for i = 1:length(Domain_list)
         in_links(i) = sum(matrix_Kd(i, :) > 0 & matrix_Kd(i, :) < Interaction_threshold);
         
         Adaptor_peptides = ismember(Peptide_list(:,1), Domain_list(i,2));
         out_links(i) = sum(sum(matrix_Kd(:, Adaptor_peptides) > 0 & matrix_Kd(:, Adaptor_peptides) < Interaction_threshold));
         peptide_count(i) = sum(Adaptor_peptides);
         normalized_out_links(i) = out_links(i) / peptide_count(i);       
end



adaptor_oncogenes = ismember(Domain_list(:,2), Oncogene_list(:,1));

[~, sort_i] = sort(in_links); [~, in_rank] = sort(sort_i);
[~, sort_i] = sort(out_links); [~, out_rank] = sort(sort_i);

isSH2 = strcmp(Domain_list(:,6), 'SH2');


p_in = ranksum(in_links(adaptor_oncogenes), in_links(~adaptor_oncogenes));
p_out = ranksum(out_links(adaptor_oncogenes), out_links(~adaptor_oncogenes));

p_sh2_in = ranksum(in_links(isSH2), in_links(~isSH2));
p_sh2_out = ranksum(out_links(isSH2), out_links(~isSH2));


figure;
plot(in_rank(isSH2), out_rank(isSH2), 'o', 'MarkerFaceColor', [.8 .8 1], 'MarkerEdgeColor', 'b', 'MarkerSize', 4); title('In Ranks vs Out Rank');
hold on;
plot(in_rank(~isSH2), out_rank(~isSH2), 'o','MarkerFaceColor', 'y', 'MarkerEdgeColor', [1 .7 0], 'MarkerSize', 4); title('In Ranks vs Out Rank');
hold off;
texobj = text(in_rank, out_rank, Domain_list(:,3));
set(texobj(adaptor_oncogenes), 'Color', [1 0 0]);
set(texobj(~adaptor_oncogenes), 'Color', [.5 .5 .5]);
axis square;
print('-depsc2', 'figure/4/AdaptorInvOut');

%Creates classification for box plot
graph_label = cell(length(Domain_list(:,1)),1);
for i = 1:length(graph_label)
    if adaptor_oncogenes(i)
        graph_label{i} = 'Oncogene';
    else
        graph_label{i} = 'Not Oncogene';
    end
end

%Plots box plots
figure; boxplot(in_rank,graph_label, 'orientation', 'horizontal'); title('Adaptor In Rank'); xlabel({[num2str(median(in_links(adaptor_oncogenes))), ' vs ', num2str(median(in_links(~adaptor_oncogenes)))] ; ['p value : ', num2str(p_in, '%4.2e')]});
set(gca,  'XTick', [], 'XLim', [0 121]);  print('-depsc2', 'figure/5/AdaptorIn');

figure; boxplot(out_rank,graph_label, 'orientation', 'horizontal'); title('Adaptor Out Rank'); xlabel({[num2str(median(out_links(adaptor_oncogenes))), ' vs ', num2str(median(out_links(~adaptor_oncogenes)))] ; ['p value : ', num2str(p_out, '%4.2e')]});
set(gca,  'XTick', [], 'XLim', [0 121]); print('-depsc2', 'figure/5/AdaptorOut');

% figure; boxplot(peptide_count,graph_label); title('Adaptor Peptide Count'); xlabel({[num2str(median(peptide_count(adaptor_oncogenes))), ' vs ', num2str(median(peptide_count(~adaptor_oncogenes)))]; ['p value : ', num2str(p_pep_count, '%4.2e')]});
% figure; boxplot(normalized_out_links,graph_label); title('Adaptor Out Connectivity Normalized for Peptides'); xlabel({[num2str(median(normalized_out_links(adaptor_oncogenes))), ' vs ', num2str(median(normalized_out_links(~adaptor_oncogenes)))]; ['p value : ', num2str(p_normal_out, '%4.2e')]});

figure; boxplot(in_rank,Domain_list(:,6), 'orientation', 'horizontal'); title('PTB vs SH2 In Rank'); xlabel({[num2str(median(in_links(isSH2))), ' vs ', num2str(median(in_links(~isSH2)))] ; ['p value : ', num2str(p_sh2_in, '%4.2e')]});
set(gca,  'XTick', [], 'XLim', [0 121]); print('-depsc2', 'figure/5/SH2PTBIn');
figure; boxplot(out_rank,Domain_list(:,6)); title('PTB vs SH2 Out Rank'); xlabel({[num2str(median(out_links(isSH2))), ' vs ', num2str(median(out_links(~isSH2)))] ; ['p value : ', num2str(p_sh2_out, '%4.2e')]});
set(gca,  'XTick', [], 'YLim', [0 121]); print('-depsc2', 'figure/5/SH2PTBOut');



%%For input into Interactive Tree of Life for fig 4a
xlswrite('/output/domain_links', [in_links in_rank out_links out_rank]); xlswrite('/output/domain_labels', Domain_list(:,3));
end

