function [ newClusterIndex ] = mergeClustersKSTest( data,clusterIndex, alpha, distanceMetric )
%mergeClustersKSTest Merges clusters based on distance distribution and KS
%test significance
% Inputs:
% x: m x n matrix, where m is the number of data points (cells) and n
% number of features
% clusterIndex: current indices of clusters
% alpha: Alpha value of significance
% distanceMetric: Metric for distance computation

if(isempty(data))
    disp('Empty matrix');
    newClusterIndex = [];
    return
end
if(nargin < 2)
    disp('Not enough arguments');
    return;
elseif(nargin < 3)
    alpha = .01;
    distanceMetric = 'euclidean';
elseif(nargin < 4)
    distanceMetric = 'euclidean';
end


maxNumIterations = 50; % Maximum number of iterations for merging

for k = 1:maxNumIterations
    if(numel(unique(clusterIndex)) <= 1)
        disp('There is only one cluster');
        break;
    end
    clusterIndex = resetIndex(clusterIndex);
    uniqueClusters = unique(clusterIndex);
    fprintf('#Iteration: %i #Clusters %d\n',k,numel(uniqueClusters));
    clsMap = zeros(numel(uniqueClusters),2);
    clsMap(:,1) = uniqueClusters;
%     fprintf('Number of Clusters %i\n',numel(uniqueClusters));
    clusterMedians = zeros(numel(uniqueClusters),size(data,2));
    for iClusters = 1:numel(uniqueClusters)
        %         ii = clusterIndex==uniqueClusters(iClusters);
        clusterMedians(iClusters,:) = median(data(clusterIndex==uniqueClusters(iClusters),:));
    end
    
    % Find nearest cluster based on median data
    [clusNN, D] = knnsearch(clusterMedians,clusterMedians,'K',2,'Distance',distanceMetric);
    clusNN = uniqueClusters(clusNN(:,2));
    clsMap(:,2) = clusNN;
    [~,ix] = sort(D(:,2),'ascend');
    clsMap = clsMap(ix,:);
    
    
    % Remove duplications
    %     maxc = size(clsMap,1);
    logicalFlag = true(size(clsMap,1),1);
    for i = 1:size(clsMap,1)
        if(~logicalFlag(i))
            continue;
        end
        ix = sum(or(clsMap(:,1:2) == clsMap(i,1),...
            clsMap(:,1:2) == clsMap(i,2)),2)>0;
        logicalFlag(ix,1)= false;
        logicalFlag(i) = true;
    end
    clsMap = clsMap(logicalFlag,:);
    m = size(clsMap,1);
    mergeClusters = true(m,1);
    threshold = alpha/m;
    for iClusters = 1:m
        numCls1= clusterIndex==clsMap(iClusters,1);
        numCls2= clusterIndex==clsMap(iClusters,2);
        
        a1 = pdist2(data(numCls1,:),clusterMedians(clsMap(iClusters,1),:),distanceMetric);
        a2 = pdist2(data(numCls1,:),clusterMedians(clsMap(iClusters,2),:),distanceMetric);
        pp1 = kstest2(a1,a2,'alpha',threshold);
%         [~,pp1] = ranksum(a1,a2,'alpha',threshold);
        a1 = pdist2(data(numCls2,:),clusterMedians(clsMap(iClusters,1),:),distanceMetric);
        a2 = pdist2(data(numCls2,:),clusterMedians(clsMap(iClusters,2),:),distanceMetric);
        pp2 = kstest2(a1,a2,'alpha',threshold);
%         [~,pp2] = ranksum(a1,a2,'alpha',threshold);
        
        if(and(pp1,pp2))
            mergeClusters(iClusters) = false;
        end
        %          [~,pvalue(iClusters)] = ranksum(a1,a2);
    end
    if(sum(mergeClusters)>0)
        for i = 1:m
            if(~mergeClusters(i))
                continue;
            end
            clusterIndex(clusterIndex==clsMap(i,1)) = clsMap(i,2);
        end
    else
        disp('No More Merging possible');
        break;
    end
    if(k==maxNumIterations)
        disp('Maximum iterations reached');
    end
end
newClusterIndex = clusterIndex;
fprintf('Number of clusters: %d\n',numel(unique(newClusterIndex)));
end

function [newIndex]  = resetIndex(indices)
uIndx = unique(indices);
newIndex = zeros(numel(indices),1);
for i = 1:numel(uIndx)
    newIndex(indices == uIndx(i)) = i;
end


end
