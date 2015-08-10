function [Links_list Linkends_sh2 Linkends_ptb] = build_circos_input_single(Adaptor, Domain_list, Peptide_list, replicate_Merged_vectors, interaction_Kd)
%Greg Koytiger June 2012
%Converts Adaptor-Adaptor interaction data into a format that can be
%rendered by Circos. Colors links by affinity and stacks them such that
%highest affinity links are drawn on top of lower affinity links

%%Clusters proteins by sequence
[Protein] = organizeProteins(Domain_list);

%%Generates "karayotype" corresponding to SH2/PTB containing proteins
Protein_length = zeros(size(Protein,1), 1); 
Protein_labels = cell(size(Protein,1),7);

for i = 1:size(Protein,1)
    Protein_length(i) = length(Domain_list{find(strcmp(Protein(i), Domain_list(:,2)), 1, 'first'),11});
    Protein_labels(i,:) = ['chr' '-' Protein(i) Protein(i) num2str(1) num2str(Protein_length(i)) 'black'];
end

%%Generates "chromosome bands" that correspond to domains
Band_labels = cell(0,7);
for i = 1:size(Protein,1)
    idx = strcmp(Protein(i), Domain_list(:,2));
    bands = cell2mat(Domain_list(idx,9:10));
    color = strrep(Domain_list(idx,6), 'SH2', 'blue');
    color = strrep(color, 'PTB', 'yellow');
    
    if(size(bands,1) == 1)
        start_band = [1 bands(1) - 1];
        end_band = [bands(2) + 1 Protein_length(i)];
        Domain_bands = [start_band; bands; end_band];
        
        temp_labels = [['band' Protein(i) 'band1', 'band1' num2cell(start_band) 'black']; 
            ['band' Protein(i) 'band2' 'band2' num2cell(bands) color]; 
            ['band' Protein(i) 'band3' 'band3' num2cell(end_band) 'black']];
        
        eliminate = (Domain_bands(:,2) - Domain_bands(:,1)) < 1;
        temp_labels(eliminate,:) = [];
        Band_labels = [Band_labels; temp_labels];
    else
        start_band = [1, min(bands(:,1)) - 1];
        mid_band = [min(bands(:,2)) + 1, (max(bands(:,1))-1)];
        end_band = [max(bands(:,2)) + 1, Protein_length(i)];
        Domain_bands = [start_band; bands; mid_band; end_band];
        temp_labels = [['band' Protein(i) 'band1' 'band1' num2cell(start_band) 'black']; 
            ['band' Protein(i) 'band2' 'band2' num2cell(bands(1,:)) color{1}];
            ['band' Protein(i) 'band3' 'band3' num2cell(mid_band) 'black'];
            ['band' Protein(i) 'band4' 'band4' num2cell(bands(2,:)) color{2}];
            ['band' Protein(i) 'band5' 'band5' num2cell(end_band) 'black']];
        eliminate = (Domain_bands(:,2) - Domain_bands(:,1)) < 0;
        temp_labels(eliminate,:) = [];
        Band_labels = [Band_labels; temp_labels]; %#ok<*AGROW>
    end
     
end
xlswrite('output/circos/sh2_karayotype.xls', [Protein_labels; Band_labels]);


%%Builds links of protein-protein interactions
replicate_Merged_vectors = replicate_Merged_vectors(ismember(replicate_Merged_vectors(:,1), cell2mat(Domain_list(:,1))) & ismember(replicate_Merged_vectors(:,2), cell2mat(Peptide_list(:,3))) & replicate_Merged_vectors(:,3) < interaction_Kd, :);

isBound = false(size(Domain_list,1),1); %%

Links_list = cell(size(replicate_Merged_vectors,1)+1000, 5); k = 1;
for i = 1:size(replicate_Merged_vectors,1)
    domain_idx = cell2mat(Domain_list(:,1)) == replicate_Merged_vectors(i,1);
    pep_idx = find(cell2mat(Peptide_list(:,3)) == replicate_Merged_vectors(i,2));
    for j = 1:length(pep_idx)
        if(strcmp(Peptide_list(pep_idx(j),1), Adaptor) || (ismember(Peptide_list(pep_idx(j),1), Domain_list(:,2)) && strcmp(Domain_list(domain_idx,2), Adaptor)))
            
            isBound(domain_idx) = true;
            Links_list{k,1} = [Peptide_list{pep_idx(j), 1}, '-', Peptide_list{pep_idx(j), 2}, '_', Domain_list{domain_idx,3}];
            Links_list{k,2} = Peptide_list{pep_idx(j), 1};
            Links_list{k,3} = str2double(Peptide_list{pep_idx(j), 2}(2:end))-5;
            Links_list{k,4} = str2double(Peptide_list{pep_idx(j), 2}(2:end))+5;
            Links_list{k,5} = ['color=(',num2str(floor(colorassign(replicate_Merged_vectors(i,3)*1000) .* [255 255 255]), '%i,%i,%i'), ')', ',z=',num2str(floor((1-replicate_Merged_vectors(i,3))*100))];
            
            Links_list{k+1,1} = [Peptide_list{pep_idx(j), 1}, '-', Peptide_list{pep_idx(j), 2}, '_', Domain_list{domain_idx,3}];
            Links_list{k+1,2} = Domain_list{domain_idx,2};
            Links_list{k+1,3} = Domain_list{domain_idx,9};
            Links_list{k+1,4} = Domain_list{domain_idx,10};
            Links_list{k+1,5} = ['color=(',num2str(floor(colorassign(replicate_Merged_vectors(i,3)*1000) .* [255 255 255]), '%i,%i,%i'), ')', ',z=',num2str(floor((1-replicate_Merged_vectors(i,3))*100))];
            
            k = k+2;
        end
    end
end

Links_list(cellfun('isempty', Links_list(:,1)),:) = [];

Linkends_sh2 = Domain_list(isBound & strcmp(Domain_list(:,6), 'SH2'), [2,9:10]);
Linkends_ptb = Domain_list(isBound & strcmp(Domain_list(:,6), 'PTB'), [2,9:10]);

xlswrite(['output/circos/', Adaptor, '_links.xls'], Links_list);  

if(~isempty(Linkends_sh2))
    xlswrite(['output/circos/', Adaptor, '_Linkends_sh2.xls'], Linkends_sh2);
end
if(~isempty(Linkends_ptb))
    xlswrite(['output/circos/', Adaptor, '_Linkends_ptb.xls'], Linkends_ptb);
end


end