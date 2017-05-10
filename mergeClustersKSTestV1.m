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
    uniqueClusters = unique(clusterIndex);
    fprintf('Number of Clusters %i\n',numel(uniqueClusters));
    clusterMedians = zeros(numel(uniqueClusters),size(data,2));
    for iClusters = 1:numel(uniqueClusters)
%         ii = clusterIndex==uniqueClusters(iClusters);
        clusterMedians(iClusters,:) = median(data(clusterIndex==uniqueClusters(iClusters),:));
    end

    % Find nearest cluster based on median data
    clusNN = knnsearch(clusterMedians,clusterMedians,'K',2,'Distance',distanceMetric);
    clusNN = clusNN(:,2);
    actualNNClsNumber = uniqueClusters(clusNN);
    pvalue  = zeros(numel(uniqueClusters),1);
    for iClusters = 1:numel(uniqueClusters)
        numCls1= sum(clusterIndex==uniqueClusters(iClusters));
        numCls2= sum(clusterIndex==actualNNClsNumber(iClusters));
        if(numCls1 <= numCls2) % Try to use the smaller cluster for distance analysis
            a1 = pdist2(data(clusterIndex==uniqueClusters(iClusters),:),clusterMedians(iClusters,:),distanceMetric);
            a2 = pdist2(data(clusterIndex==uniqueClusters(iClusters),:),clusterMedians(clusNN(iClusters),:),distanceMetric);
        else
            a1 = pdist2(data(clusterIndex==actualNNClsNumber(iClusters),:),clusterMedians(iClusters,:),distanceMetric);
            a2 = pdist2(data(clusterIndex==actualNNClsNumber(iClusters),:),clusterMedians(clusNN(iClusters),:),distanceMetric);
        end
        [~,pvalue(iClusters)] = kstest2(a1,a2,alpha);
%          [~,pvalue(iClusters)] = ranksum(a1,a2);
    end
    threshold = alpha/(numel(uniqueClusters)); % Threshold with bonferroni correction
    if(sum(pvalue > threshold) > 0)        
        cls1 = uniqueClusters(pvalue > threshold);
        cls2 = actualNNClsNumber(pvalue > threshold);
        uNumCls = unique(cls2);
        for i = 1:numel(uNumCls)
            l2 = cls1(cls2 == uNumCls(i));
            for j = 1:numel(l2)
                if(numel(unique(clusterIndex)) > 1)
                    fprintf('Merging clusters %i and %i\n',l2(j),uNumCls(i));
                    clusterIndex(clusterIndex==(l2(j))) = uNumCls(i);
                end
            end
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

end

