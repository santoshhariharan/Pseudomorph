function [ newIndex ] = assignDataToNearestClusters( data,indexToBeAssigned,clusterIndex )
%assignDataToNearestClusters Assigns data to nearest clusters

newIndex = zeros(size(data,1),1);
newIndex(~indexToBeAssigned) = clusterIndex;


% dataAlreadyAssigned = data(~indexToBeAssigned,:);
uClusterIndex = unique(clusterIndex);
clusterCenters = zeros(numel(uClusterIndex),size(data,2));
for iCls = 1:numel(uClusterIndex)
    clusterCenters(iCls,:) = mean(data(newIndex==uClusterIndex(iCls),:)); 
end
idx = knnsearch(clusterCenters,data(indexToBeAssigned,:));
newIndex(indexToBeAssigned) = uClusterIndex(idx);
end
