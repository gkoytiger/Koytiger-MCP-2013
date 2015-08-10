%This script performs all the fitting and data analysis for protein
%microarrays
%Greg Koytiger June 2011
%Alexis Kaushansky, Ethan Karp  (basis of Reformat code) May 2009

close all 

%load('Data/Unfit_Data.mat');  
load('Data/Domain_list.mat');
load('Data/Peptide_list.mat');
load('Data/sanger_oncogenes.mat');
load('Data/bushman_oncogenes.mat');
load('Data/uniprot_RTK.mat')
 
%Load to skip titration fitting (Comment out if using Unfit_Data)
load('Data/Merged_data.mat');

%%%%%%%%%%%These are the variables that can be changed %%%%%%%%%%%%%%%%%%
useParallel = true;         %Enables parallelization of code (takes ~16hrs on quad core in parallel mode to fit all titrations)
PercentExclude = 25;        %threshold for excluding outlier data using trimmean function

threshold_Kd_low = 10^-2;   %lower bound of valid in vitro Kd (uM) (lowest titration point)
threshold_Kd_high = 5;      %upper bound of valid in vitro Kd (uM) (top titration point)
threshold_R2 = 0.90;        %Lowest acceptable R^2
threshold_FOB = 0.90;       %Lowest acceptable fold over background (expressed as a percentile)
 
interaction_Kd = 1;         %Threshold for 'in vivo' interactions (uM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fits microarray images to generate titrations, requires Unfit_Data

generate_figure_1a(Domain_list, PercentExclude, interaction_Kd, Data);

%%Fits saturation binding curves to microarray images
Fit_Data = cell(size(Data,1), 1);

if(useParallel)
    matlabpool open; %Enable parallel computation
end

for PlateRow = 1:321
   MaxConcentrations = Data{PlateRow,2};
   ArrayType = Data{PlateRow,4};
   
   RawData_list = Reformat(Data, PlateRow, ArrayType, PercentExclude);  
   [AverageData_list Cy5NormalizedData] = DataAverage(RawData_list, PercentExclude); 
   Fit_Data{PlateRow} = Fit(Cy5NormalizedData, AverageData_list, MaxConcentrations, PercentExclude);
end

clear PercentExclude PlateRowstart PlateRowend PlateRow RawData_list AverageData_list MaxConcentrations Cy5NormalizedData ArrayType

if(useParallel)
    matlabpool close; %closes parallel threads
end
 
%%Merges all data and filters into format:
%%Domain_id Peptide_id Kd R2 FoldOverBackground Kd_95CI_low Kd_95CI_high

Merged_data = vertcat(Fit_Data{:});
save('Data/Merged_data.mat', 'Merged_data', '-v7.3');

clear Data Fit_Data

%%Filters data based on Kd, R2, Signal to Noise
threshold_Merged_vectors = threshold_Merged_data(Merged_data, threshold_Kd_high, threshold_Kd_low, threshold_R2, threshold_FOB);


%%Merge data replicates
replicate_Merged_vectors = replicate_Merged_data(threshold_Merged_vectors);
 
matrix_Kd = generate_matrix_Kd(replicate_Merged_vectors, Peptide_list, Domain_list);

%%Analyzes relationship between RTK connectivity and oncogencity (visualized in figure 2a)
[RTK_p RTK_conn_onc RTK_conn_not_onc sorted_RTKs] = generateFigure2a(matrix_Kd, Peptide_list, RTK_list, Oncogene_list, Interaction_threshold);

%%Generates connectivity diagram of RTKs (basis for Fig 1B, 1C, 2B
generateConnectivityDiagrams( matrix_Kd, Peptide_list, Domain_list, sorted_RTKs, interaction_Kd )


%%Generate input for Circos to visualize adaptor binding (basis for Fig 3a,3b, 3d)
build_circos_input(Domain_list, Peptide_list, replicate_Merged_vectors, interaction_Kd);
build_circos_input_single('GRB2',Domain_list, Peptide_list, replicate_Merged_vectors, interaction_Kd);
build_circos_input_single('PTPN11',Domain_list, Peptide_list, replicate_Merged_vectors, interaction_Kd);
build_circos_input_single('SHC1',Domain_list, Peptide_list, replicate_Merged_vectors, interaction_Kd);

%%Generates the rank diagram for adaptors (Fig 4b)
generateFigure4b(matrix_Kd, Peptide_list, Domain_list, bushman_oncogenes, interaction_Kd);

%%Generate Ouput tables
headers = {'Protein', 'pY Site', 'Peptide Sequence', 'Gene Id', 'Uniprot Id'};
xlswrite('Output/Supplemenatry_Table_S1.xls', [headers; Peptide_list(:, [1 2 4 5 6])]);
 
headers ={'Domain', 'Gene Id', 'Uniprot Id', 'Domain Type', 'Cloned Start', 'Cloned End'};
xlswrite('Output/Supplemenatry_Table_S2.xls', [headers; Domain_list(:, 3:8)]);

Interaction_list = build_output(replicate_Merged_vectors, Peptide_list, Domain_list, interaction_Kd);
headers = {'Adaptor', 'Gene Id', 'Uniprot Id', 'Protein', 'pY Site', 'Gene Id', 'Uniprot Id', 'Peptide Sequence', 'Kd (uM)'};
xlswrite('Output/Supplemenatry_Table_S3.xls', [headers; Interaction_list]);
  

%%Estimates false positive rate (fpr) and false discovery rate (fdr) by assuming that the least connected domains (bottom 20%) are inactive
%%in our assay and that all binding events for these domains are stochastic false positives
neg_percentile= .20;
[ fpr fdr estimate_false_hits] = estimate_false_positives( Interaction_list, Domain_list, Peptide_list, neg_percentile, Merged_data );

clear useParallel PercentExclude threshold_Kd_low threshold_Kd_high threshold_R2 threshold_FOB  headers 
