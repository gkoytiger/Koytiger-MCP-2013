function [matrix_Kd] = generate_matrix_Kd(Merged_vectors, Peptides, Domain_info)
%Builds a matrix of the kd values where matrix_Kd(Domain,Peptide) = Kd in
%nM

%Domain_Numbers = cell2mat(Domain_info(:,1));
Domain_Numbers = cell2mat(Domain_info(:,1));
Peptide_Numbers = cell2mat(Peptides(:,3));

matrix_Kd = zeros(size(Domain_info,1),size(Peptides,1));

for Index = 1:size(Merged_vectors,1)
    matrix_Kd((Merged_vectors(Index,1) == Domain_Numbers), (Merged_vectors(Index,2) == Peptide_Numbers)) = (Merged_vectors(Index,3)*1000);
end


end

