% Plot clusters with PCA:
% For each cluster, use 10 percent of cluster population
% closest to cluster median

% Use bubble plot minsize(5), masize(10) of the data to represent medians
% Use kruskal's algorithm to plot MST of the medians
clc;close all;
% cntrl = {'ergic';'golgin';'tabik';'CB5-ER';'CYTO';'ER-PRO';'LAMINA';'MAO';...
%             'METAXIN';'MITO-CCO';'PEROXISOME-1';'PQC-PSS1';'RAB5A';'RAB7A';...
%             'VAMP1A'};
% cntrl = {'ergic';'golgin';'tabik';'ER-PRO'};
% uniqueClusterIndex = unique(mergingIndex(:,iIterations+1));
sampleDataforViewing =  false;
samplingNumber = 2000;
%% Remove clusters With 5% or less of any control
minControlPercentage = .05;

uniqueClusterIndex = unique(mergingIndex(:,iIterations+1));
cDistribution = zeros(numel(uniqueClusterIndex),numel(cntrl));
numberCells = cell2mat(clusterFileMap(:,2));
for i = 1:numel(uniqueClusterIndex)
    ii = mergingIndex(:,iIterations+1) == uniqueClusterIndex(i);
    
    for j = 1:numel(cntrl)
        jj = strcmpi(clusterFileMap(:,3),cntrl{j,1});
        cDistribution(i,j) = sum(numberCells(and(ii,jj),1));
    end
end

clustersToBeRemoved = bsxfun(@rdivide,cDistribution,sum(cDistribution));
clustersToBeRemoved = clustersToBeRemoved<minControlPercentage;
clustersToBeRemoved = sum(clustersToBeRemoved,2)==numel(cntrl);
%% Consider only contorls supplied
controlsToKeep = false(size(clusterFileMap,1),1);
for i = 1:numel(cntrl)
    ii = strcmpi(clusterFileMap(:,3),cntrl{i,:});
    controlsToKeep(ii) = true;
end
clusterFileMap = clusterFileMap(controlsToKeep,:);
% Get points near cluster centroids

dataForPCA = [];%Bad programming - change later
numData = 0;
clusterCentroid = zeros(numel(uniqueClusterIndex),sum(finalFeatures));
numPointsPerCluster = zeros(numel(uniqueClusterIndex),1);
% Start iterations over clusters to get pooled data
for i = 1:numel(uniqueClusterIndex)
%     Get files for every cluster
    ii = mergingIndex(:,iIterations+1) == uniqueClusterIndex(i);
    fNames = clusterFileMap(ii,:);
    fprintf('In cluster: %d\n',uniqueClusterIndex(i));
    d=[];
    for j = 1:size(fNames,1)
        fprintf('-- %s\n',fNames{j,1});
        fid = fopen(fNames{j,1},'r');
        data = fread(fid,[fNames{j,2} sum(finalFeatures)],'double');
        fclose(fid);
        d = [d;data];
    end
%     clusterCentroid(i,:) = median(d);
    clusterCentroid(i,:) = mean(d);
%     idx = randperm(size(d,1),50);
    if(sampleDataforViewing)
        idx = knnsearch(d,clusterCentroid(i,:),'k',2000);
        dataForPCA = [dataForPCA;d(idx',:)];
        numData = numData + numel(idx);
    else
        dataForPCA = [dataForPCA;d];
    end
    numPointsPerCluster(i) = size(d,1);
    fprintf('--Original size: %d Samsize: %d\n',size(d,1),round(.05.*size(d,1)));
    
end
dataForPCA = [clusterCentroid;dataForPCA];
dataForPCA = clusterCentroid;
disp('done');
clear  numData i ii fNames
clear d j fid data idx
%% Use principal components
numdims = 2;
numRedType = 'PCA';
% [coeff] = pca(dataForPCA,'Centered',true);
clusterCentroidPCA = compute_mapping(dataForPCA,numRedType,numdims);
clusterCentroidPCA = clusterCentroidPCA(1:numel(uniqueClusterIndex),:);

clusterCentroidPCA = clusterCentroidPCA(~clustersToBeRemoved,:);
numPointsPerCluster = numPointsPerCluster(~clustersToBeRemoved);
uniqueClusterIndex = uniqueClusterIndex(~clustersToBeRemoved);
clusterCentroid = clusterCentroid(~clustersToBeRemoved,:);
cDistribution = cDistribution(~clustersToBeRemoved,:);
% Plot the cluster centroid with size of the marker proportional to size
% percentage of data points
pp = round(sqrt(numPointsPerCluster));
pp = pp - min(pp);
sl = (5./(max(pp)-min(pp)));

% pp(pp~=1) = pp(pp~=1)+eps;
% pp(pp==1) = 1 - eps;
markerSize = round((pp.*sl) + 4);
figure;hold on;
for i = 1:size(clusterCentroidPCA,1)
%     mrkerSize
    if(numdims==3)
        plot3(clusterCentroidPCA(i,1),clusterCentroidPCA(i,2),clusterCentroidPCA(i,3),'or',...
            'Markerfacecolor',[.6 .6 .6],'MarkerSize',markerSize(i,1),...
            'Markeredgecolor','none','LineStyle','None');
    else
        plot(clusterCentroidPCA(i,1),clusterCentroidPCA(i,2),'or',...
            'Markerfacecolor',[.6 .6 .6],'MarkerSize',markerSize(i,1),...
            'Markeredgecolor','none','LineStyle','None');
    end
end
hold off;
% xlabel('t-SNE-1');ylabel('t-SNE-2');
title([numRedType ' Cluster centroid']); 

%% Describe data
numberCellsPerControl = zeros(size(cntrl,1),1);
numberClustersPerControl = zeros(size(cntrl,1),1);
for j = 1:numel(cntrl)
    jj = strcmpi(clusterFileMap(:,3),cntrl{j,1});
    numberCellsPerControl(j,1) = sum(numberCells(jj,1));
    numberClustersPerControl(j,1) = sum(jj);
end
% Plot number of cells per cluster
% figure;
% semilogy(1:size(clusterCentroidPCA,1),sum(cDistribution,2),'or','MarkerFaceColor','r');
% hold on;
% semilogy(1:size(clusterCentroidPCA,1),repmat(200,...
%             1,size(clusterCentroidPCA,1)),'-b');
% % semilogy(1:numel(uIndx),repmat(400,1,numel(uIndx)),'-k');
% hold off;
% ylabel('Cells per cluster');xlabel('Cluster number');
%% Perform minimum spanning tree 
bkgColor = [.4 .4 .4];
am = ones(numel(uniqueClusterIndex));
for i = 1:numel(uniqueClusterIndex)
    am(i,i) = 0;
end
D = pdist2(clusterCentroidPCA,clusterCentroidPCA,'cityblock');
% idx = knnsearch(clusterCentroidPCA,clusterCentroidPCA,'k',5);
% D = knn2jaccard(idx(:,2:end));
% D = full(D);
% D = D+D';
[w,xst] = kruskal(am,D);

figure;hold on;grid on;
for i = 1:size(clusterCentroidPCA,1)
%     mrkerSize
if(numdims==3)
    plot3(clusterCentroidPCA(i,1),clusterCentroidPCA(i,2),clusterCentroidPCA(i,3),'or',...
        'Markerfacecolor',bkgColor,'MarkerSize',markerSize(i,1),...
        'Markeredgecolor','none','LineStyle','None'); 
%     text(clusterCentroidPCA(i,1),clusterCentroidPCA(i,2),...
%         clusterCentroidPCA(i,3),num2str(i))
else
    plot(clusterCentroidPCA(i,1),clusterCentroidPCA(i,2),'or',...
        'Markerfacecolor',bkgColor,'MarkerSize',markerSize(i,1),...
        'Markeredgecolor','none','LineStyle','None'); 
%     text(clusterCentroidPCA(i,1),clusterCentroidPCA(i,2),num2str(i))
    
end
end

for i = 1:size(xst,1)    
    coord = clusterCentroidPCA(xst(i,:),:)';
    if(numdims==3)
        line(coord(1,:),coord(2,:),coord(3,:),'Linewidth',.2,'Color',[.4 .4 .4]);
    else
        line(coord(1,:),coord(2,:),'Linewidth',.2,'Color',[.4 .4 .4]);
    end
end
hold off;
% xlabel('t-SNE-1');ylabel('t-SNE-2');
% title('Cluster centroid'); 
title([numRedType ' Cluster centroid']); 
%%
viewMSTPie2(clusterCentroidPCA,cDistribution,jet(3),cntrl,xst);


%% Plot cluster distribution & enrichment

% Get cluster distrbution with controls



% pval = enrich_score(cDistribution);
% pval = pval < (.05./(numel(uniqueClusterIndex).*numel(cntrl)));

% For each control highlight enriched clusters in red
% for j = 1:numel(cntrl)
%     
%     figure;hold on;
%     for i = 1:size(clusterCentroidPCA,1)
%         %     mrkerSize
%         if(pval(i,j))
%             if(numdims==3)
%                 plot3(clusterCentroidPCA(i,1),clusterCentroidPCA(i,2),clusterCentroidPCA(i,3),'or',...
%                 'Markerfacecolor',[1 0 0],'MarkerSize',markerSize(i,1),...
%                 'Markeredgecolor','none','LineStyle','None');
%                 
%             else
%                 plot(clusterCentroidPCA(i,1),clusterCentroidPCA(i,2),'or',...
%                 'Markerfacecolor',[1 0 0],'MarkerSize',markerSize(i,1),...
%                 'Markeredgecolor','none','LineStyle','None');
%             end
%         else
%             if(numdims==3)
%                 plot3(clusterCentroidPCA(i,1),clusterCentroidPCA(i,2),clusterCentroidPCA(i,3),'or',...
%                 'Markerfacecolor',[.8 .8 .8],'MarkerSize',markerSize(i,1),...
%                 'Markeredgecolor','none','LineStyle','None');
%             else
%                 plot(clusterCentroidPCA(i,1),clusterCentroidPCA(i,2),'or',...
%                 'Markerfacecolor',[.8 .8 .8],'MarkerSize',markerSize(i,1),...
%                 'Markeredgecolor','none','LineStyle','None');
%             end
%         end
%         
%         
%     end
%     hold off;
%     xlabel('t-SNE-1');ylabel('t-SNE-2');
%     title([numRedType ' ' cntrl{j,:}]); 
%     title([numRedType ' Cluster centroid']); 
    %     pval(:,j)
   
% end




