function [Proteins] = organizeProteins(Domain_list)
%Greg Koytiger July 2012
%Organizes domains around the outside of rectangular visualization of
%binding.

%Identifies all unique proteins
%Note: SH2 has hierarchy over PTB domains, thus if protein has both PTB and
%SH2 domain (eg SHC) its SH2 domain location takes precedence

[~, m1, ~] = unique(Domain_list(:,2), 'first');
Protein_list = Domain_list(m1,:);
[~,idx] = sort(m1);
Protein_list = Protein_list(idx,:);
j=1;k=1;sh2_sequences = cell(1,1); ptb_sequences = cell(1,1);
sh2_labels = cell(1,1); ptb_labels = cell(1,1);

for i = 1:size(Protein_list,1)

    if(strcmp('SH2', Protein_list{i,6}))
        sh2_sequences{j} = Protein_list{i,11}(Protein_list{i,9}: Protein_list{i,10});
        sh2_labels{j} = Protein_list{i,2};
        j = j+1;
    else
        ptb_sequences{k} = Protein_list{i,11}(Protein_list{i,9}: Protein_list{i,10});
        ptb_labels{k} = Protein_list{i,2};
        k = k+1;
    end
    
end

sh2Tree = seqlinkage(seqpdist(sh2_sequences), 'average', sh2_labels);
ptbTree = seqlinkage(seqpdist(ptb_sequences), 'average', ptb_labels);



%Retrieves the order of the sequence linkage
Proteins = [get(sh2Tree, 'LeafNames'); get(ptbTree, 'LeafNames')];

end