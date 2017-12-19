function [ h ] = viewMSTPieCluster2Cluster(x,clusterIndex,treatmentIndex,...
                                    map,adjacencyMatrix,edgeWidth,...
                                    edgeColor,treatmenNames)
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
edgeWidthVal = (v - min(v))./(max(v)-min(v));
% v = edgeWidthVal;
edgeWidthVal = 1-edgeWidthVal;
edgeWidthVal = 4*edgeWidthVal+.001;
xst = [r c];



% In the same figure plot data using patches - Uses inbuilt Matlab function
% logic
% 
% numSmallPatch = 30;
numTreatments = size(disMat,2);
numSteps = 30;
h=figure;hold all;

for i = 1:size(xst,1)
    coord = x(xst(i,:),:)';  
    cl = edgeColor(xst(i,1),xst(i,2),:);
    line(coord(1,:),coord(2,:),'Linewidth',edgeWidthVal(i),'Color',cl);    
end

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
    text(cCenter(1,1),cCenter(1,2),num2str(iClusters),'Color',[0 0 0]);
%     legend(uYhat)
end
hold off;axis fill
% legend(hall,treatmentNames);
% clear hall
axis off;
end
