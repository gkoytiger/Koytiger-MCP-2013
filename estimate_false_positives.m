function [ fpr fdr total_false_hits] = estimate_false_positives( Interaction_list, Domain_list, Peptide_list, neg_percentile, Merged_data )
%Greg Koytiger June 2012
%Assumes bottom X percentile of domains are inactive and that any binding
%is stochastic false positives (fpr). Uses this to estimate how much of the data
%is liable to be false positives (fdr). As this measure is fairly sensitive
%to percentile threshold, it should be considered an estimate of true fpr
%and fdr

%Finds all titrations performed on peptides and domains of interest
filtered_data = Merged_data(ismember(Merged_data(:,1), cell2mat(Domain_list(:,1))) & ismember(Merged_data(:,2), cell2mat(Peptide_list(:,3))), 1);

domain_hits = zeros(size(Domain_list,1), 1);
total_hits = size(Interaction_list,1);

%Tallies the total number of hits for each domain
for i = 1:size(Domain_list,1)
    domain_hits(i) = sum(strcmp(Domain_list(i,3), Interaction_list(:,1)));
end

[~, idx] = sort(domain_hits); [~, idx] = sort(idx); %Transforms hits into ranks
neg_idx = ceil(neg_percentile * size(Domain_list,1)); %Transforms percentile cut off into rank cut off

false_hits = sum(domain_hits(idx < neg_idx)); %Extracts total hits corresponding to domains below rank cut off
total_titrations = sum(ismember(filtered_data, cell2mat(Domain_list(idx < neg_idx,1))) );   %Finds the total number of titrations performed on negative domains

fpr = false_hits / total_titrations; 
total_false_hits = fpr * size(filtered_data,1);
fdr = total_false_hits / total_hits;

end

