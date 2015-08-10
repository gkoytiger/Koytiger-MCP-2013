function [Background_list] = findBackground(RawData_list, Wells, PercentExclude)
%Finds the appropriate background for each spot
Background_Spot_Even = 500; 
Background_Spot_Odd = 501;


Background_Prefetch_Even = RawData_list((RawData_list(:,1) == Background_Spot_Even), :);
Background_Prefetch_Odd  = RawData_list((RawData_list(:,1) == Background_Spot_Odd), :);


Background_list = zeros(length(RawData_list),1);

parfor listI = 1:size(RawData_list,1)
    if mod(Wells(listI),2) == 0 %determine which background spot to use if even or odd well
       Background_list(listI) = Spot_Average(Background_Spot_Even,RawData_list(listI,2),RawData_list(listI,3), Background_Prefetch_Even, PercentExclude);

    else
       Background_list(listI) = Spot_Average(Background_Spot_Odd,RawData_list(listI,2),RawData_list(listI,3), Background_Prefetch_Odd, PercentExclude);

    end
end