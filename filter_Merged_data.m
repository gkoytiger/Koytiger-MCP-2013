function [filtered_Merged_vectors] = filter_Merged_data(filtered_Merged_vectors, Peptide_list, Domain_list)
%Greg Koytiger June 2011

Domains = cell2mat(Domain_list(:,1));
Peptides = cell2mat(Peptide_list(:,3));

eliminate = ~ismember(filtered_Merged_vectors(:,1), Domains) | ~ismember(filtered_Merged_vectors(:,2), Peptides);
filtered_Merged_vectors(eliminate,:) = [];
 
end


