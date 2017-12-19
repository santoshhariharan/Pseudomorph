function [ h ] = viewMSTPie2(x,clusterIndex,treatmentIndex,...
                                    map,adjacencyMatrix,edgeWidth,...
                                    treatmentNames)
%viewMSTPie2 Plots data in x using grouping in y with colours of yhat
% Plots cluster centroids based on x. 
% Input:
% x - Reduced dimensional plot
% clusterIndex - Cluster Indices
% treatmentIndex - Treatment Indices
% map - Color map for treatments
% adjacencyMatrix - Matrix of MST connection
% edgeWidth - Width of the edges
% xst - Minimum Spanning tree

if(size(x,2)>2)
    error('@viewScatterPie: Only two dmensional plot allowed');    
end

% Compute distribution matrix
uniqueClusterIndex = unique(clusterIndex);
uniqueTreatmentIndex = unique(treatmentIndex);
uniqueClusterIndex = uniqueClusterIndex(uniqueClusterIndex~=0);
disMat = zeros(numel(uniqueClusterIndex),numel(uniqueTreatmentIndex));
for iCls = 1:numel(uniqueClusterIndex)
    ii = clusterIndex == uniqueClusterIndex(iCls);
    for jTreatment = 1:numel(uniqueTreatmentIndex)
        jj = treatmentIndex == uniqueTreatmentIndex(jTreatment);
        disMat(iCls,jTreatment) = sum(ii.*jj);
    end
end

% Estimate radius for each cluster
axisRange = abs(max(x)-min(x));
maxRadius = .04*min(axisRange);
minRadius = .01*min(axisRange);
rsize = sum(disMat,2);
rsize = (rsize - min(rsize))./(max(rsize)- min(rsize));
rsize = (maxRadius - minRadius)*rsize + minRadius;

% Get Edges
adjacencyMatrix = tril(adjacencyMatrix);
[r,c] = find(adjacencyMatrix);
lv = find(adjacencyMatrix);
[v] = edgeWidth(lv);
if(numel(v)>1)
edgeWidthVal = (v - min(v))./(max(v)-min(v));
v = edgeWidthVal;
edgeWidthVal = 1*edgeWidthVal+.1;
else
    edgeWidthVal = 1;
end
xst = [r c];



% In the same figure plot data using patches - Uses inbuilt Matlab function
% logic
% 
% numSmallPatch = 30;
numTreatments = size(disMat,2);
% angleStart = 0;
numSteps = 30;
% rsize = 5;
% hold on;
h=figure;hold all;
% 1-(v(i)*
% Plot points
points2Plot = unique([r;c]);

for i = 1:size(xst,1)
    coord = x(xst(i,:),:)';    
    line(coord(1,:),coord(2,:),'Linewidth',edgeWidthVal(i),'Color',[.7 .7 .7]);    
end
% for i = 1:points2Plot
%     cl = map(treatmentIndex(points2Plot(i)),:);
    scatter(x(points2Plot,1),x(points2Plot,2),45,...
        map(treatmentIndex(points2Plot),:),'filled'); 
% end




for iClusters = 1:size(disMat,1)
    startAngle = 0;  
    treatProportion = disMat(iClusters,:)./sum(disMat(iClusters,:));
    cCenter = x(uniqueClusterIndex(iClusters),:);
    hall=[];
    for jTreatment =  1:numTreatments
        step = treatProportion(jTreatment)*2*pi/numSteps;
        theta = [0 0:step:treatProportion(jTreatment)*2*pi 0];
        theta = theta + startAngle;
        [x1, y1] = pol2cart(theta,[0 rsize(iClusters).*(ones(1,numel(theta)-2)) 0]);
        x1 = x1 + cCenter(1,1);
        y1 = y1 + cCenter(1,2);
        h1 = patch(x1,y1,jTreatment,'FaceColor',map(jTreatment,:),'EdgeColor','none');
        hall = [hall;h1];
        startAngle = max(theta);        
    end    
%     legend(uYhat)
end

%     scatter(x(points2Plot,1),x(points2Plot,2),45,...
%         map(treatmentIndex(points2Plot),:),'filled','MarkerEdgeColor',[0 0 0]);
hold off;axis tight
% legend(hall,treatmentNames);
% clear hall
axis off;
end

