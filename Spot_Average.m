function [SpotAverage] = Spot_Average(Domain_number, Peptide_number, Peptide_concentration, Data_list, PercentExclude)
%Returns the trimmed mean average of a spot 

%working_averages = zeros(1,4);
%Indices = intersect(intersect(find(Data_list(:,1) == Domain_number),find(Data_list(:,2)==Peptide_number)),find(Data_list(:,3)==Peptide_concentration));

working_averages = Data_list(Data_list(:,1) == Domain_number & Data_list(:,2) == Peptide_number & Data_list(:,3) == Peptide_concentration, 4); %this is a 4x3 matrix with the quadruplicates down the columns


%Takes the trimmed mean of the replicates
SpotAverage = trimmean(working_averages,PercentExclude);


            
    

