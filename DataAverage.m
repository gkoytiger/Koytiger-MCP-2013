function [AverageData_list Cy5NormalizedData] = DataAverage(RawData_list, PercentExclude)
%Greg Koytiger June 2011
%Transforms RawData by subtracting background and finds spot averages for
%all quadruplicates

AverageData_list = zeros(length(unique(RawData_list(:,1:3),'rows')),5);
AverageData_list(:,1:3) = unique(RawData_list(:,1:3),'rows'); %this finds all of the unique interacting pairs

%Removes blanks and errors with peptides = 9999 

BadPeptideIndex = AverageData_list(:,2) == 9999;
AverageData_list(BadPeptideIndex,:) = [];

BackgroundSubtractedData = RawData_list(:,1:6);
Cy5NormalizedData = RawData_list(:,1:6);

%Finds Signal to Noise of individual spots
SpotSignalToNoise = [RawData_list(:,1:3) (RawData_list(:,4) ./ RawData_list(:,6))];


%Background Subtraction
BackgroundSubtractedData(:, 4) = RawData_list(:,4) - RawData_list(:,6);


%Data Cy5 Normalize
Cy5NormalizedData(RawData_list(:,5) == 0, :) = [];
Cy5NormalizedData(:, 4) = BackgroundSubtractedData(:,4) ./ RawData_list(:,5);

%Uncomment to remove background subtraction
%Cy5NormalizedData(:, 4) = Cy5NormalizedData(:,4) ./ RawData_list(:,5);



%Average Cy3/Cy5 Values and average signal to noise of spots
for listI = 1:length(AverageData_list) 
   AverageData_list(listI,4) = Spot_Average(AverageData_list(listI,1), AverageData_list(listI,2), AverageData_list(listI,3), Cy5NormalizedData, PercentExclude);
   AverageData_list(listI,5) = Spot_Average(AverageData_list(listI,1), AverageData_list(listI,2), AverageData_list(listI,3), SpotSignalToNoise, PercentExclude);
end

end