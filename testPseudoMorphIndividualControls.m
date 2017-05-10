% Load sampled data
clear;
clc;
pth='F:\Projects\Proteinlocalization\PseudoMorph\TestDataForPseudoMorph';
load(fullfile(pth,'sampleDataCellFeat.mat'));
load(fullfile(pth,'parameters.mat'));
% Perform feature reduction once with PCA based and once with unsupervised
% greedy optimization
[redFeatures, featureScores ]= unsupervisedPCASelect(data,30);
data=data(:,redFeatures);
% % Community Detection
% % Perform community detection over multiple values of K
% % Use Jaccard distance & then euclidean distance
% numNeighbours = [5:5:300];
% numCls = zeros(numel(numNeighbours),1);
% for jNeighbors = 1:numel(numNeighbours)
%     indx = phenograph( data, numNeighbours(jNeighbors));
%     numCls(jNeighbors) = numel(unique(indx));
% end
% clc;
% addpath(genpath('F:\Projects\Proteinlocalization\PseudoMorph\Code\drtoolbox'));
% score = compute_mapping(data, 't-SNE', 2);
% score = score(:,1:2);
% rmpath(genpath('F:\Projects\Proteinlocalization\PseudoMorph\Code\drtoolbox'));
%% Go over every control

allIndex = zeros(size(data,1),1);
uC = unique(textdata);
setK = 50;cnt = 0;
for iCtrl = 1:numel(uC)
    ii = strcmpi(textdata,uC{iCtrl,:});
    nData = data(ii,:);
    indx = phenograph( nData, setK);
    allIndex(ii) = indx+cnt;
    cnt = cnt+numel(unique(indx))+1;
    fprintf('Iteration%i\n',iCtrl);
end

% Compute cluster centroids & Iteratively merge clusters
maxClsIter = 50;
uIndex = unique(allIndex);
clsCent = zeros(numel(uIndex),size(data,2));
mergeTree = zeros(numel(uIndex),maxClsIter);
for iCls = 1:numel(uIndex)
    ii = allIndex == uIndex(iCls);
    clsCent(iCls,:) = mean(data(ii,:));
end

%% Search Nearest Clusters by Centroids
maxMergingIter = 50;
mergingIndex = zeros(numel(uIndex),maxMergingIter+1);
mergingIndex(:,1) = [1:numel(uIndex)];
for iIterations = 1:maxClsIter
    if(iIterations ~= 1)        
        nClsCent = zeros(numel(uIndex),size(data,2));
        for i = 1:numel()
    else
        nClsCent = clsCent;
    end
    [idx, dist] = knnsearch(clsCent,clsCent,'K',2);
%     idx = idx(:,2);
%     dist = dist(:,2);
    
    logicalFlag = true(size(idx,1),1);
    for i = 1:size(idx,1)
        if(~logicalFlag(i))
            continue;
        end
        ix = sum(or(idx(:,1:2) == idx(i,1),...
            idx(:,1:2) == idx(i,2)),2)>0;
        logicalFlag(ix,1)= false;
        logicalFlag(i) = true;
    end
    mergeIdx = idx(logicalFlag,:);
    
    for i =1:size(mergeIdx,1)
        ix1Data = data(allIndex==mergeIdx(i,1),:);
        ix2Data = data(allIndex==mergeIdx(i,2),:);
%         perfromPermutationMerging
        c1 = median(ix1Data);
        c2 = median(ix2Data);
        pd1 = pdist2(ix1Data,c1);
        pd2 = pdist2(ix2Data,c1);
        pp1 = ranksum(pd1,pd2,'alpha',.05);
        pd1 = pdist2(ix1Data,c2);
        pd2 = pdist2(ix2Data,c2);
        pp2 = ranksum(pd1,pd2,'alpha',.05);
        mergeVal = ~and(pp1,pp2);
        if(mergeVal)
            fprintf('--Merging clusters %i and %i\n',idx(i,1),idx(i,2));
            ii = mergingIndex(:,iIterations) == idx(i,1);
            mergingIndex(ii,iIterations+1) = idx(i,2);
        end        
    end
%    allIndex =  resetIndexValues(allIndex,mergingIndex(:,iIterations),...
%                         mergingIndex(ii,iIterations+1));
end




% Start Merging Process
% disp('Ready');