% Try & merge clusters formed previously using KS test
%

load(fullfile(pwd,'ClusterOutput','clusterParameters.mat'),'clusterFileMap','finalFeatures');
% Read individual files and compute cluster medians
% cntrl = {'ergic';'golgin';'tabik';'CB5-ER';'CYTO';'ER-PRO';'LAMINA';'MAO';...
%             'METAXIN';'MITO-CCO';'PEROXISOME-1';'PQC-PSS1';'RAB5A';'RAB7A';...
%             'VAMP1A'};
% cntrl = {'golgin';'ergic';'ER-PRO'};
for i = 1:1000
    if(i>size(clusterFileMap,1))
        break;
    end
    if(isempty(clusterFileMap{i,3}))
        break;
    end
end
clusterFileMap=clusterFileMap(1:i-1,:);
%% Consider only contorls supplied
controlsToKeep = false(size(clusterFileMap,1),1);
for i = 1:numel(cntrl)
    ii = strcmpi(clusterFileMap(:,3),cntrl{i,:});
    controlsToKeep(ii) = true;
end
clusterFileMap = clusterFileMap(controlsToKeep,:);
%%
% clusterCentroid = zeros(size(clusterFileMap,1),sum(finalFeatures));
% for i = 1:size(clusterFileMap,1)
%     fid = fopen(clusterFileMap{i,1},'r');
%     data = fread(fid,[clusterFileMap{i,2} sum(finalFeatures)],'double');
%     fclose(fid);
%     clusterCentroid(i,:) = median(data);
% end
% idx = knnsearch(clusterCentroid,clusterCentroid,'K',2,'Distance','Euclidean');
% idx = idx(:,2);

% Number of iterations
clc;
maxMergingIter = 50;
mergingIndex = zeros(size(clusterFileMap,1),maxMergingIter+1);
mergingIndex(:,1) = 1:size(clusterFileMap,1);

for iIterations = 1:maxMergingIter
    fprintf('Iteration: %i\n',iIterations);
    uniqueClusters = unique(mergingIndex(:,iIterations));
    clusterCentroid = zeros(numel(uniqueClusters),sum(finalFeatures));
    for jUniqueClusters = 1:numel(uniqueClusters)
        ii = mergingIndex(:,iIterations) == uniqueClusters(jUniqueClusters);
        fNames = clusterFileMap(ii,:);
%         numRows = sum(cell2mat(clusterFileMap(ii,2)))
        d  = [];
        for k = 1:size(fNames,1)
            fid = fopen(fNames{k,1},'r');
            data = fread(fid,[fNames{k,2} sum(finalFeatures)],'double');
            fclose(fid);
            d = [d;data];
        end
        clusterCentroid(jUniqueClusters,:) = median(d);
%         clusterCentroid(jUniqueClusters,:) = mean(d);
    end
    
%     Check for closest cluster for each one
    [idx,d] = knnsearch(clusterCentroid,clusterCentroid,'K',2,'Distance','Euclidean');
%     idx = idx(:,2);
    d = d(:,2);
%     Remove multiple assignment for any cluster
    [~,ii] = sort(d);
    idx = idx(ii,:);
    clear d;
    
    
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
    idx = idx(logicalFlag,:);
    mergeVal = true(size(idx,1),1);
    alpha = .01./size(idx,1);
    mergingIndex(:,iIterations+1) = mergingIndex(:,iIterations);
    for i = 1:size(idx,1)
%         disp(i)
        i1 = mergingIndex(:,iIterations) == idx(i,1);
        i2 = mergingIndex(:,iIterations) == idx(i,2);
        fNames = clusterFileMap(i1,1:2);
        mergeVal(i) = testClusterMerging( clusterFileMap(i1,1:2),clusterFileMap(i2,1:2),...
            sum(finalFeatures),alpha );
        if(mergeVal(i))
            fprintf('--Merging clusters %i and %i\n',idx(i,1),idx(i,2));
            ii = mergingIndex(:,iIterations) == idx(i,1);
            mergingIndex(ii,iIterations+1) = idx(i,2);
        end
        
    end
    mergingIndex(:,iIterations+1) = reindexClusters(mergingIndex(:,iIterations+1));
    if(sum(mergeVal)==0)       
%         fprintf('Final number of clusters: %d\n',numel(unique(mergingIndex(:,iIterations+1))));
        disp('No more Merging possible');
        fprintf('Final number of clusters: %d\n',numel(unique(mergingIndex(:,iIterations+1))));
        break;
        
    end
    
%     disp('Round1');

    
    
end
clear data alpha fid fNames i i1 i2 idx ii
clear ix jUniqueClusters k logicalFlag maxMergingIter mergeVal uniqueClusters


