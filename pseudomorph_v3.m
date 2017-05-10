% Define Global paramaters
clc;clear;
maxSplits = 5;
numSplits = maxSplits;
numNeighbours = 21;
maxPointsPerSplit = 10000;
distanceMetric = 'Euclidean';
maxSampleSize = 100000;
computeWellMedianFeatures = false;
% clc;
load('parameters.mat');
x = getappdata(0,'alldata');
fid = fopen(fullfile(param.rootpath,'_controls_temp.bin'),'r');
x.data = fread(fid,[size(x.textdata,1) sum(param.datafeat)],'double');
fclose(fid);
eleNames = {'Golgin';'Ergic';'TABIK'};
groups = getGroupIndices(x.textdata(:,9),eleNames);
uGrps = unique(groups);
indexKeep = false(size(x.data,1),1);
featureIndex = zeros(sum(param.datafeat),1);
fprintf('Input paramaters:\n Num-Neighbors: %i\n Distance Metric: %s\n Maximum points per split: %i\n',numNeighbours,distanceMetric,maxPointsPerSplit);
%%

if(computeWellMedianFeatures)
    % Compute Well Medians and get features for well medians
    x.wellMedianData = zeros(800,sum(param.datafeat));
    x.wellMedianTextData = cell(800,size(x.textdata,2));
    wellColumn = 1;
    plateColumn = 2;
    uniquePlates = unique(x.textdata(:,plateColumn));
    uniqueWells = unique(x.textdata(:,wellColumn));
    count = 1;
    disp('Computing well medians');
    % fprintf('Number of Completed wells completed \t');
    fprintf('Number of Completed wells completed ... %4d%%',0);
    for iGroups = 1:numel(uGrps)
        ii = groups == uGrps(iGroups);
        for jPlates = 1:numel(uniquePlates)
            jj = strcmpi(x.textdata(:,plateColumn),uniquePlates{jPlates,:});
            for kWells = 1:numel(uniqueWells)
                kk = strcmpi(x.textdata(:,wellColumn),uniqueWells{kWells,:});
                if(sum(ii.*jj.*kk)>0)
                    x.wellMedianData(count,:) = median(x.data(logical(ii.*jj.*kk),param.datafeat));
                    tmp = x.textdata(logical(ii.*jj.*kk),:);
                    x.wellMedianTextData(count,:) = tmp(1,:);
                    fprintf('\b\b\b\b\b%4d%%',count);
                    count = count+1;
                end
            end
        end
    end
    fprintf('\n');
    indexToRemove = sum(x.wellMedianData==0,2) >= sum(param.datafeat);
    x.wellMedianData = x.wellMedianData(~indexToRemove,:);
    x.wellMedianTextData = x.wellMedianTextData(~indexToRemove,:);
    wellGroups = getGroupIndices(x.wellMedianTextData(:,9),eleNames);
    f = featureReduction(x.wellMedianData,wellGroups);
    if(~isempty(f))
        featureIndex(f) = 1;
    end
end

% Sample from x to keep maxSampleSize from each group
cellsToKeep = zeros(size(x.data,1),1);
for iGroups = 1:numel(uGrps)
    ii = find(groups == uGrps(iGroups));
    cellsToKeep(ii(randperm(numel(ii),min(numel(ii),maxSampleSize)))) = 1;
end 
cellsToKeep = logical(cellsToKeep);
% Keep only randomly sampled data
% ydata = x.data(cellsToKeep,param.datafeat);
ydata = x.data(cellsToKeep,:);
ytextdata = x.textdata(cellsToKeep,:);

clear x cellsToKeep

% Clear variables

clear wellColumn plateColumn uniquePlates uniqueWells count tmp ii jj kk iGroups jPlates kWells indexToRemove x.data x.textdata
%% Split data into 5 Parts & compute Clustering separately


[m,n] = size(ydata);
groups = getGroupIndices(ytextdata(:,9),eleNames); 
newGroups = ones(m,1);
% Random splitting into groups
pp = randperm(m);
newGroups(pp(1:floor(m/2))) = 2;
clear pp;
% newGroups = groups;
allDataIndex = zeros(m,1);
numClusters = 0;
uniqueGroups = unique(newGroups);
numNeighbours = [21:10:51];
maxPointsPerSplit = [10000 25000];
alphaVal = [.05 .01 .001 .0001 .00001];
numClusterKS = zeros(1,numel(alphaVal));
for jNeighbours = 1:numel(numNeighbours)
    for kPointsperSplit = 1:numel(maxPointsPerSplit)
        for iGroups = 1:numel(uniqueGroups)
            fprintf('Running group %i\n',iGroups);
            groupIndex = newGroups == uniqueGroups(iGroups); % Get only groups for sampled data
            if(sum(groupIndex) > maxPointsPerSplit(kPointsperSplit)) % Split data if the number of data points per groups is greater than a set value
                %         splitIndex = false(sum(groupIndex),1);
                numSplits = round(sum(groupIndex)/maxPointsPerSplit(kPointsperSplit));
                cpar = cvpartition(sum(groupIndex),'kfold',numSplits);
                splitIndex = false(numel(groupIndex),numSplits);
                tmp = find(groupIndex);
                for iSplits = 1:numSplits
                    splitIndex(tmp(cpar.test(iSplits)),iSplits) = true;
                end
            else
                numSplits = 1;
                splitIndex = true(sum(groupIndex),1);
            end
            fprintf('Size of group %i\n',sum(groupIndex));
            
            for iSplits = 1:numSplits
                fprintf('  In split %i\n',iSplits);
                tStart = tic;
                idx = knnsearch(ydata(logical(splitIndex(:,iSplits).*groupIndex),:),...
                    ydata(logical(splitIndex(:,iSplits).*groupIndex),:),'K',numNeighbours(jNeighbours),'Distance',distanceMetric);
                idx = idx(:,2:end);
                sim = knn2jaccard(idx);
                [com] = cluster_jl(sim,1,0);
                [~,maxMod] = max(com.MOD);
                indx = com.COM{1,maxMod}';
                indx = indx + numClusters;
                numClusters = numClusters + numel(unique(indx));
                allDataIndex(and(splitIndex(:,iSplits),groupIndex),1) = indx;
                fprintf('    Time in seconds %.1f\n',toc(tStart));
                fprintf('    Found %i clusters in split %i and groups %i\n',numel(unique(indx)),iSplits,iGroups);
            end
            
        end
        %
        for i = 1:numel(alphaVal)
            allDataIndexCpy = allDataIndex;
            valClusters = allDataIndex ~= -1;
            allDataIndexCpy(valClusters) = mergeClustersKSTest(ydata(valClusters,:),allDataIndex(valClusters),alphaVal(i),distanceMetric);
            fprintf('#Neighbors %d #PointsPerSplit %d #BeforeMerging %d #AfterMerging %d AlphaVal %f \n',numNeighbours(jNeighbours),maxPointsPerSplit(kPointsperSplit),...
                                                                                                    numel(unique(allDataIndex(allDataIndex~=-1))),...
                                                                                                    numel(unique(allDataIndexCpy(allDataIndex~=-1))));    
            allDataIndexCpy(valClusters) = removeOutlierClusters(ydata(valClusters,:),allDataIndexCpy(valClusters),distanceMetric,.001,.001);
            valClusters = allDataIndexCpy ~= -1;
            numClusterKS(i) = numel(unique(allDataIndexCpy(valClusters)));
        end
        titleStr = sprintf('KS Merging: #Neighbors %d #PointsPerSplit %d',numNeighbours(jNeighbours),maxPointsPerSplit(kPointsperSplit));
        figure;semilogx(alphaVal,numClusterKS,'-xr');title(titleStr);xlabel('Alpha value');ylabel('Number of clusters');
    end
end
return;
disp('Done');

%% Check distribution
% uniqueGroups = unique(groups);
uniqueGroups = unique(newGroups);
uniqueClusters = unique(allDataIndexCpy);
clusterDistribution = zeros(numel(uniqueClusters),numel(uniqueGroups));
for iClusters = 1:numel(uniqueClusters)
    if(uniqueClusters(iClusters) == -1)
        continue;
    end
    ii = allDataIndexCpy == uniqueClusters(iClusters);
    for jGroups = 1:numel(uniqueGroups)
        jj = groups == uniqueGroups(jGroups);
        clusterDistribution(iClusters,jGroups) = sum(ii.*jj); 
    end
end
return;
%% Cluster distribution based on plate
uniqueGroups = unique(ytextdata(:,2));
uniqueClusters = unique(allDataIndex);
clusterDistributionPlate = zeros(numel(uniqueClusters),numel(uniqueGroups));
for iClusters = 1:numel(uniqueClusters)
    if(uniqueClusters(iClusters) == -1)
        continue;
    end
    ii = allDataIndex == uniqueClusters(iClusters);
    for jGroups = 1:numel(uniqueGroups)
        jj = strcmpi(ytextdata(:,2),uniqueGroups{jGroups,:});
        clusterDistributionPlate(iClusters,jGroups) = sum(ii.*jj); 
    end
end
disp('Done')

return;
%%
% Sample data from each cluster and plot a t-SNE 
% grpSplitData = []; grpSplitIndex = [];grpClsIndex = [];
% c = cvpartition(groups,'holdout',.05);
% Cluster medians
uniqueClusters = unique(allDataIndex);
clusterMedians = zeros(numel((uniqueClusters)),sum(param.datafeat));
for i = 1:numel(uniqueClusters)
    clusterMedians(i,:) = median(x.resampledData(allDataIndex == uniqueClusters(i),:)) ;
end
index = [];
for i = 1: numel(uGrps)
    ii = groups == uGrps(i);
    for j = 1:numel(uniqueClusters)
        jj = allDataIndex == uniqueClusters(j);
        tmp = find(ii.*jj);        
        index = [index ;tmp(randperm(numel(tmp),min(numel(tmp),100)))];
    end
end

ggrps = groups(index);
cclsIndex = allDataIndex(index);    
% y = compute_mapping(clusterMedians, 'Sammon', 2);
% [~,y] = princomp(clusterMedians);
[coeff,y] = princomp(x.resampledData(index,:));
cMediansReduced = clusterMedians*coeff(:,[1 2]);
y = y(:,1:3);
% [~,y] = princomp(
% Perform SDA
% map = colormap(hsv(20));close all;
% marker1 = ['s','d','+'];
% figure; hold on;
% for i = 1: size(y,1)
%     plot(y(i,1),y(i,2),'Marker',marker1(1),'MarkerFaceColor',map(i,:),'LineStyle','none',...
%                         'MarkerSize',5,'MarkerEdgeColor','none');
%                     legendNames = [legendNames;cellstr(['Cluster-' num2str(i)])];
% end
% hold off;
map1 = [1 0 0; 0 1 0; 0 0 1];
uggrps = unique(ggrps);
figure; hold on;
for i = 1:numel(uggrps)
    ii = ggrps == uggrps(i);
    plot3(y(ii,1),y(ii,2),y(ii,3),'Marker','o','MarkerFaceColor',map1(i,:),'LineStyle','none',...
                        'MarkerSize',3,'MarkerEdgeColor','none')
end
hold off;
%% f = featureReduction(x.data(allIdxCpy>0,param.datafeat),allIdxCpy(allIdxCpy>0));grpData);
%
% map1 = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; .5 .5 .5; .75 .2 .3];
map = colormap(hsv(20));close all;
marker1 = ['d','s','+'];
% [~,maxMod] = max(com.MOD);
% indx = com.COM{1,maxMod}';
% uIndx = unique(indx);
uggrps = unique(ggrps);
ucclsIndex = unique(cclsIndex);
figure; hold on;
legendNames={};
count = 1;
for i = 1: numel(uggrps)
    ii = ggrps == uggrps(i);
    for j = 1: numel(ucclsIndex)
        jj = cclsIndex == ucclsIndex(j);
        if(sum(ii.*jj) > 0)
            plot(y(logical(ii.*jj),1),y(logical(ii.*jj),2),'Marker',marker1(i),'MarkerFaceColor',map(count,:),'LineStyle','none',...
                        'MarkerSize',5,'MarkerEdgeColor','none');
    %     plot3(y(indx==i,1),y(indx==i,2),y(indx==i,3),'Marker','o','MarkerFaceColor',map(i,:),...
    %                         'LineStyle','none','MarkerSize',5,'MarkerEdgeColor','none');
                        legendNames = [legendNames;cellstr(['Cluster-' num2str(j)])];
                        count = count +1;
                    
        end
         
    end
end

for i = 1:numel(ucclsIndex)
    plot(cMediansReduced(i,1),cMediansReduced(i,2),'Marker','o','MarkerFaceColor',map(i,:),'LineStyle','none',...
                        'MarkerSize',7,'MarkerEdgeColor','none');
end
legend(legendNames,'Interpreter','None');
hold off;
%%

A = ones(size(cMediansReduced,1),size(cMediansReduced,1));
for i = 1:size(A,1)
    A(i,i) = 0;
end
zz = pdist2(cMediansReduced(),cMediansReduced,'cosine');
[C,D,F] = kruskal(A,zz);
legendNames=[];
figure;hold on
for i = 1:size(cMediansReduced,1)
    plot(cMediansReduced(i,1),cMediansReduced(i,2),'Marker','o','MarkerFaceColor',map(i,:),'LineStyle','none',...
                        'MarkerSize',8,'MarkerEdgeColor','none');
                    legendNames = [legendNames;cellstr(['Cluster-' num2str(i)])];
end
for i = 1:size(D,1)
    line([cMediansReduced(D(i,1),1) cMediansReduced(D(i,2),1)],[cMediansReduced(D(i,1),2) cMediansReduced(D(i,2),2)],'LineWidth',2,'Color',[.8 .8 .8]);
%     pause;
end
legend(legendNames,'Interpreter','None');
hold off
%% Save each

