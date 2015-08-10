%% Greg Koytiger Aug 2012
%% Requires CGDS Matlab Toolbox

cgdsURL = 'http://www.cbioportal.org/public-portal/';
cancerStudies = getcancerstudies(cgdsURL);


geneticProfilesData = cell(5, length(cancerStudies.cancerTypeId));

geneticProfilesLabel = {'rna_seq_mrna'};


for i = 1:length(cancerStudies.cancerTypeId)
    
    geneticProfiles = getgeneticprofiles(cgdsURL, cancerStudies.cancerTypeId{i});
    
    %Collects genetic profiling data of interest
    for k = 1:length(geneticProfilesLabel)
       
        loc = find(strcmp(geneticProfiles.geneticProfileId, [cancerStudies.cancerTypeId{i},'_',geneticProfilesLabel{k}]));
        if(~isempty(loc))
          profileData = getprofiledata(cgdsURL, strcat(cancerStudies.cancerTypeId{i},'_all'),geneticProfiles.geneticProfileId{loc}, Proteins(:,1), true);
          geneticProfilesData{k,i} = profileData.data;
        end
    end
    

end
publicData = {'coadread_tcga', 'kirc_tcga', 'brca_tcga'};
cancers = cancerStudies.cancerTypeId;

geneticProfilesData(~ismember(cancers, publicData),:) = [];
mergedData = [geneticProfilesData{:}];
mergedData(:, isnan(mergedData(1,:))) = [];
mergedData(isnan(mergedData(:,1)),:) = [];

%Figure Sup 2A
figure; hist(sum(mergedData > 1,1)/143,45);
%Figure Sup 2B
figure; hist(sum(mergedData > 1,2)/1437,45);

