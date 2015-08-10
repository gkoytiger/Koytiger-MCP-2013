function Interaction_list = build_output(replicate_Merged_vectors, Peptide_list, Domain_list, interaction_Kd)
%Greg Koytiger May 2012
%Builds a table of all the measured interactions

%Filters out all data that does not correspond to peptides and domains of
%interest
replicate_Merged_vectors = replicate_Merged_vectors(ismember(replicate_Merged_vectors(:,1), cell2mat(Domain_list(:,1))) & ismember(replicate_Merged_vectors(:,2), cell2mat(Peptide_list(:,3))) & replicate_Merged_vectors(:,3) < interaction_Kd, :);

Interaction_list = cell(size(replicate_Merged_vectors,1)+1000, 8);
k = 1;

for i = 1:size(replicate_Merged_vectors,1)
    domain_idx = cell2mat(Domain_list(:,1)) == replicate_Merged_vectors(i,1);
    pep_idx = find(cell2mat(Peptide_list(:,3)) == replicate_Merged_vectors(i,2));
 
    for j = 1:length(pep_idx)
        Interaction_list{k,1} = Domain_list{domain_idx,3};
        Interaction_list{k,2} = Domain_list{domain_idx,4};
        Interaction_list{k,3} = Domain_list{domain_idx,5};
    
        Interaction_list{k,4} = Peptide_list{pep_idx(j), 1};
        Interaction_list{k,5} = Peptide_list{pep_idx(j), 2};
        Interaction_list{k,6} = Peptide_list{pep_idx(j), 5};
        Interaction_list{k,7} = Peptide_list{pep_idx(j), 6};
        Interaction_list{k,8} = Peptide_list{pep_idx(j), 4};
    
        Interaction_list{k,9} = replicate_Merged_vectors(i, 3);
        
        k = k+1;
    end
end

Interaction_list(cellfun('isempty', Interaction_list(:,1)),:) = [];

end

