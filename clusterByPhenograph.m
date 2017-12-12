function [ clusterIndex,clusterCentroid,centroidFraction, uniqueIndex] = clusterByPhenograph( x,k )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

indx = phenograph(x,k,'graphtype','Jaccard');
uIndx = unique(indx);
clusterCentroid = zeros(numel(uIndx),size(x,2));
centroidFraction = zeros(numel(uIndx),1);
uniqueIndex= zeros(numel(uIndx),1);
clusterIndex = zeros(numel(indx),1);
for j = 1:numel(uIndx)
    jj = find(indx == uIndx(j));
    idx = knnsearch(x(jj,:),mean(x(jj,:)),'K',1);
    idx = idx(end);
    idx = jj(idx);
    clusterCentroid(j,:) = x(idx,:);
    uniqueIndex(j,1) = idx;
    centroidFraction(j,1) = numel(jj)/numel(indx);
    clusterIndex(jj) = idx;
end
uniqueIndex = sort(uniqueIndex);
end

