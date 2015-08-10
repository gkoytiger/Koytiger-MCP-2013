function generateConnectivityDiagrams( matrix_Kd, Peptide_list, Domain_list, sorted_RTKs, interaction_Kd )
%Greg Koytiger May 2012
%Generates rectangular diagram for visualizing RTK connectivity
pgap = 100;     %Spacing between peptides
diameter = 10;  %Circle diameter of each domain
width = 38;     %Domains across width of plot
height = 20;    %Domains across length of plot
gap = 50;       %Spacing between domains
text_offset = 40;

Domain_layout = zeros(size(Domain_list,1),2);

i = 1; countw = 0;     counth = height;

%Re-sorts Domain list by multiple sequence alignment
[Proteins] = organizeProteins(Domain_list);
[~,loc] = ismember(Domain_list(:,2), Proteins);
[~,idx] = sort(loc);
Domain_list = Domain_list(idx,:);
matrix_Kd = matrix_Kd(idx,:);

%Builds the first layer of domains in the figure
secondLayer = ~cellfun('isempty', strfind(Domain_list(:,3), '-PTB')) | ~cellfun('isempty', strfind(Domain_list(:,3), '-NC'));


while ~(countw == width)
    if(~secondLayer(i))
        
        Domain_layout(i,1) = gap * countw;
        Domain_layout(i,2) = gap * counth;
        countw = countw + 1;
    end
    
    i = i+1;
end
countw = width;

while ~(counth == 0)
    if(~secondLayer(i))
        Domain_layout(i,1) = gap * countw;
        Domain_layout(i,2) = gap * counth;
        counth = counth - 1;
    end
    
    i = i+1;
end
counth = 0;

while ~(countw == 0)
    if(~secondLayer(i))
        Domain_layout(i,1) = gap * countw;
        Domain_layout(i,2) = gap * counth;
        countw = countw - 1;
    end
    i = i+1;
end
countw = 0;

while ~(counth == height)
    if(~secondLayer(i))
        Domain_layout(i,1) = gap * countw;
        Domain_layout(i,2) = gap * counth;
        counth = counth + 1;
    end
    i = i+1;
end


%%Builds the coordinates of the second layer of domains for the figure

idx = find(secondLayer);
firstLayout = Domain_layout(~secondLayer,:);
firstLayerDomains = Domain_list(~secondLayer,:);

for i = 1:sum(secondLayer)
    Domain_layout(idx(i),:) = mean(firstLayout(strcmp(Domain_list(idx(i),2), firstLayerDomains(:,2)), :) ,1);
end

isTop = Domain_layout(:,2) == height*gap;
isBottom = Domain_layout(:,2) == 0;
isRight = Domain_layout(:,1) == width*gap & ~isTop & ~isBottom;
isLeft = Domain_layout(:,1) == 0 & ~isTop & ~isBottom;

Domain_layout(isTop & secondLayer,2) = Domain_layout(isTop & secondLayer,2) + gap;
Domain_layout(isBottom & secondLayer,2) = Domain_layout(isBottom & secondLayer,2) - gap;
Domain_layout(isRight & secondLayer,1) = Domain_layout(isRight & secondLayer,1) + gap;
Domain_layout(isLeft & secondLayer,1) = Domain_layout(isLeft & secondLayer,1) - gap;


%Creates the lines for connecting domains on same protein
Connect = zeros(4,118);
idx = 1;

for i = 1:size(Domain_layout,1)
    tf = find(strcmp(Domain_list(i,2), Domain_list(:,2)));
    
    if(length(tf) == 2)
        Connect(:, idx) = [Domain_layout(tf(1),1) Domain_layout(tf(2),1) Domain_layout(tf(1),2) Domain_layout(tf(2),2)]';
        idx = idx + 1;
    end
    
    if(size(tf,1) == 3)
        Connect(:,idx)    = [Domain_layout(tf(1),1) Domain_layout(tf(2),1) Domain_layout(tf(1),2) Domain_layout(tf(2),2)]';
        Connect(:,idx+1)  = [Domain_layout(tf(1),1) Domain_layout(tf(3),1) Domain_layout(tf(1),2) Domain_layout(tf(3),2)]';
        Connect(:,idx+2)  = [Domain_layout(tf(2),1) Domain_layout(tf(3),1) Domain_layout(tf(2),2) Domain_layout(tf(3),2)]';
        idx = idx+3;
    end
end
Connect = unique(Connect', 'rows')';


%Hides redundant text
hidetext = ismember(Domain_list(:,2), Domain_list(secondLayer,2));
hidetext(secondLayer) = false;

text_layout = Domain_layout(~hidetext,:);
labels = Domain_list(~hidetext,2);

isTop = text_layout(:,2) >= height*gap;
isBottom = text_layout(:,2) <= 0;
isRight = text_layout(:,1) >= width*gap & ~isTop & ~isBottom;
isLeft = text_layout(:,1) <= 0 & ~isTop & ~isBottom;

text_layout(isTop,2) = text_layout(isTop,2) + text_offset;
text_layout(isBottom,2) = text_layout(isBottom,2) - text_offset;
text_layout(isRight,1) = text_layout(isRight,1) + text_offset;
text_layout(isLeft,1) = text_layout(isLeft,1) - text_offset;

%%Generates peptide layout and connections

isSH2 = strcmp(Domain_list(:,6), 'SH2');


for i = 1:size(sorted_RTKs,1)
    
    figure;
    %subplot(8,5,i);
    hold on;
    
    isRTK = ismember(Peptide_list(:,1), sorted_RTKs(i,1));
    RTK_peps = Peptide_list(isRTK,:); RTK_kd = matrix_Kd(:,isRTK);
    Pep_layout = zeros(size(RTK_peps,1),2);
    Pep_layout(:,2) = height / 2 * gap; %Puts peptides in center of image
    wcenter = (width / 2 * gap);
    
    for j = 1:size(Pep_layout,1)
        Pep_layout(j, 1) = wcenter - (size(Pep_layout,1)/2-j) * pgap; %Aligns peptides around horizontal center
        binds = RTK_kd(:,j) > 0 & RTK_kd(:,j) < interaction_Kd*1000;
        aff = RTK_kd(binds,j);
        Y = Domain_layout(binds,:);
        X = repmat(Pep_layout(j,:), size(Y,1), 1);
        if ~isempty(aff)
            
            for k = 1:size(aff,1)
                line([X(k,1) Y(k,1)]',[X(k,2) Y(k,2)]', 'Color', colorassign(aff(k)), 'LineWidth', 2);
            end
            
        end
    end
    line(Connect(1:2,:), Connect(3:4,:), 'Color', 'black', 'LineWidth', 1);
    
    plot(Pep_layout(:,1),Pep_layout(:,2), 'o', 'MarkerSize', diameter, 'MarkerFaceColor', [1 .8 .8], 'MarkerEdgeColor', 'r', 'LineWidth', 2);
    
    text(Pep_layout(:,1),Pep_layout(:,2)+text_offset, RTK_peps(:,2),'HorizontalAlignment', 'left', 'FontSize', 12, 'FontWeight', 'bold', 'Rotation', 90, 'FontName', 'Arial');
    text(wcenter, height*gap+250,sorted_RTKs(i,1), 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
    axis off;
    
    plot(Domain_layout(isSH2,1),Domain_layout(isSH2,2), 'o', 'MarkerSize', diameter, 'MarkerFaceColor', [.8 .8 1], 'MarkerEdgeColor', 'b', 'LineWidth', 2);
    plot(Domain_layout(~isSH2,1),Domain_layout(~isSH2,2), 'o', 'MarkerSize', diameter, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', [1 .7 0], 'LineWidth', 2);
    
    textobj = text(text_layout(:,1), text_layout(:,2), labels);
    set(textobj(isTop) , 'HorizontalAlignment', 'left', 'FontSize', 10, 'Rotation', 90, 'FontName', 'Arial');
    set(textobj(isBottom) , 'HorizontalAlignment', 'left', 'FontSize', 10, 'Rotation', -90, 'FontName', 'Arial' );
    set(textobj(isRight) , 'HorizontalAlignment', 'left', 'FontSize', 10, 'FontName', 'Arial');
    set(textobj(isLeft) , 'HorizontalAlignment', 'right', 'FontSize', 10, 'FontName', 'Arial');
    set(gcf, 'Position', [150 150 1200 650], 'PaperPositionMode','auto');
    print('-depsc2', ['figure/Supplementary/', [num2str(i), '-', sorted_RTKs{i,1}]]);
    
    plot(Pep_layout(:,1),Pep_layout(:,2), 'o',  'MarkerSize', 1, 'MarkerFaceColor', [1 .8 .8], 'MarkerEdgeColor', 'r');
    plot(Domain_layout(isSH2,1),Domain_layout(isSH2,2),  'o', 'MarkerSize', 1, 'MarkerFaceColor', [.8 .8 1], 'MarkerEdgeColor', 'b');
    plot(Domain_layout(~isSH2,1),Domain_layout(~isSH2,2),  'o', 'MarkerSize', 1, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', [1 .7 0]);
    
    
end


end

