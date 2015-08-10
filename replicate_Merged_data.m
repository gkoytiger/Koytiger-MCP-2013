function [replicate_Merged_vectors] = replicate_Merged_data(threshold_Merged_vectors)
%Greg Koytiger June 2011

%Takes the weighted average of the log Kd values derived from replicate probings
%Assigns weights to be the multiple of the R2 times the fold over background percentile

SNRs = threshold_Merged_vectors(:,5);
R2s  = threshold_Merged_vectors(:,4);

log_Kd = log10(threshold_Merged_vectors(:,3));

unique_interactors = unique(threshold_Merged_vectors(:,1:2),'rows');

replicate_Merged_vectors = [unique_interactors zeros(length(unique_interactors),3)];


for Interaction = 1:length(replicate_Merged_vectors)
    Indexes = intersect(find(threshold_Merged_vectors(:,1) == replicate_Merged_vectors(Interaction,1)),find(threshold_Merged_vectors(:,2) == replicate_Merged_vectors(Interaction,2)));
    Kds_toaverage = log_Kd(Indexes);
    SNRs_toaverage = SNRs(Indexes);
    R2s_toaverage = R2s(Indexes);
    
    weights = SNRs_toaverage .* R2s_toaverage;
    
    replicate_Merged_vectors(Interaction,3) = exp(sum(Kds_toaverage .* weights)./ sum(weights));
    replicate_Merged_vectors(Interaction,4) = sum(R2s_toaverage);
    replicate_Merged_vectors(Interaction,5) = sum(SNRs_toaverage);   
end

end
