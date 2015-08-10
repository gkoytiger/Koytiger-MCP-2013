function [threshold_Merged_vectors] = threshold_Merged_data(threshold_Merged_vectors, threshold_Kd_high, threshold_Kd_low,threshold_R2, threshold_SignalToNoise)
%Greg Koytiger June 2011

%Excludes data that does not meet the Kd R2 and Signal to Noise thresholds
%initializes the filters where 1 equals passed filter, 0 does not pass
%filter

filters = zeros(length(threshold_Merged_vectors),3) + 1;

%Filters out fold over background over a certain percentile
threshold_Merged_vectors(:,5) = generate_percentiles(threshold_Merged_vectors(:,5));
filters(:,1) = threshold_Merged_vectors(:,5) > threshold_SignalToNoise;

%Filters based on Kd
filters(:,2) = threshold_Merged_vectors(:,3) < threshold_Kd_high;
filters(:,3) = threshold_Merged_vectors(:,3) > threshold_Kd_low;

%Filters based on R2
filters(:,4) = threshold_Merged_vectors(:,4) > threshold_R2;


not_passed_filters = sum(filters,2) < 4;
threshold_Merged_vectors(not_passed_filters, :) = [];

end


