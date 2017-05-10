% Pseudomorph: Per control
% Run pseudo morph on individual controls. 
% Run pseudomorph on the collected centroids again to group them to get the
% final set of centroids
% Assign all points to individual centroids and keep the propoertions
% algorithm steps
% Read data - Load Individual control files (cleaned data)
% Reduce the number of features to meaningful 30 -
% Sample based on local density
% Use k = 5 and create a sparse jaccard graph for phenograph
% Save the centroids (Labels need not be stored)
% Recluster the centroids using phenograph with an optimal value of k?
% For all data, assign data to nearest centroids
% Visualize with PCA (Based on sample)
%********** All data points at once
% Module - 1:
% Define Variables 
clear;clc;close all;
warning off;
numDimsVis = 2;
fprintf('Starting pseudomorph\n');
pth='F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data';
load(fullfile(pth,'parameters.mat'));% Load parameter file
param.rootpath = pth;
intensityFeature = 'Ch2_INT_Cell_intensity';
intFeat = strcmpi(intensityFeature,param.datahdr);
roiIntFeat = strcmpi('Ch2_MOR_cell_ROI_AvgIntensity',param.datahdr);
nucAreaFeat = strcmpi('Ch1_MOR_Nucleus_area',param.datahdr);
cellAreaFeat = strcmpi('Ch1_MOR_Cell_area',param.datahdr);
filePrefix = '.txt';
randPrc = .1;
% numFeatures = 30;
fNames = dir(pth);
columnForControls = 9;
columnForOrganelle = 10;
% featureReduction = true;
clustersPerLandmark = true;

maxMinTType = false;
% Module 3: Read & Load data after filtering
fprintf('Module 3.......\n');
mxRw = 1000000;
allD = zeros(mxRw,sum(param.datafeat));
allInten = zeros(mxRw,1);
allMorRatio = zeros(mxRw,1);
allTxt = cell(mxRw,1);
allTxtOrg = cell(mxRw,1);
allMorIntensity = zeros(mxRw,1);
cnt = 1;
fprintf('Completed Reading................');
for iFiles = 3:size(fNames,1)
    fprintf('\b\b\b\b\b\b\b\b\b%8.3f%%',iFiles*100./size(fNames,1));
    if(fNames(iFiles).isdir)
        continue;
    end
    tok = regexpi(fNames(iFiles).name,filePrefix,'match');
    if(isempty(tok))
        continue;
    else
        D = readfiles(cellstr(fNames(iFiles).name),param);
    end
    
    %     Remove cells out of focus
    focus = getFocusFilteredData(D.data,param);
    D.data = D.data(focus,:);
    D.textdata = D.textdata(focus,:);
    ii = (D.data(:,strcmpi('Ch1_INT_Nucleus_intensity',param.datahdr))./...
        D.data(:,strcmpi('Ch1_INT_Cytoplasm_intensity',param.datahdr)))>3.5;
    jj = (D.data(:,strcmpi('Ch1_INT_Nucleus_intensity_stddev',param.datahdr))./...
        D.data(:,strcmpi('Ch1_INT_Cytoplasm_intensity_stddev',param.datahdr)))>3.5;
    ii = and(ii,jj);
%     ii = true(size(D.data,1),1);
    D.data = D.data(ii,:);
    D.textdata = D.textdata(ii,:);
    allInten(cnt:cnt+size(D.data,1)-1,:) = D.data(:,intFeat); 
    allMorIntensity(cnt:cnt+size(D.data,1)-1,:) = D.data(:,roiIntFeat);
    allMorRatio(cnt:cnt+size(D.data,1)-1,:) = D.data(:,nucAreaFeat)./...
                                            D.data(:,cellAreaFeat);    
    allD(cnt:cnt+size(D.data,1)-1,:) = D.data(:,param.datafeat);
    allTxt(cnt:cnt+size(D.data,1)-1,:)= D.textdata(:,columnForControls);
    allTxtOrg(cnt:cnt+size(D.data,1)-1,:)= D.textdata(:,columnForOrganelle);
    cnt = cnt+size(D.data,1);
end
if(cnt<mxRw)
    allD = allD(1:cnt-1,:);
    allInten = allInten(1:cnt-1,:);
    allTxt = allTxt(1:cnt-1,:);
    allMorRatio = allMorRatio(1:cnt-1,:);
    allMorIntensity = allMorIntensity(1:cnt-1,:);
    allTxtOrg = allTxtOrg(1:cnt-1,:);
end

ii = sum(isnan(allD),2) ==0;
allD = allD(ii,:);
allInten = allInten(ii,:);
allTxt = allTxt(ii,:);
allMorRatio = allMorRatio(ii,:);
allMorIntensity = allMorIntensity(ii,:);
allTxtOrg= allTxtOrg(ii,:);
fprintf('\n');

% Low expression
ii = allInten>100 & allMorRatio <= .7;
allD = allD(ii,:);
allTxt = allTxt(ii,:);
allTxtOrg = allTxtOrg(ii,:);
allInten = allInten(ii,:);
allMorRatio = allMorRatio(ii,:);
allMorIntensity = allMorIntensity(ii,:);
fprintf('# Cells removed by intensity %i\n',sum(~ii));
% allInten = allInten(ii,:);
if(maxMinTType)
    allD = bsxfun(@minus,allD, min(allD));
    allD = bsxfun(@rdivide,allD,max(allD) - min(allD));
    allInten = (allInten - min(allInten))./(max(allInten) - min(allInten));
else
    allInten = zscore(allInten);
    allD = zscore(allD);
end
% Remove intensity correlated features
rho = corr(allD,allInten);
newHeader = param.datahdr(1,param.datafeat);
% Retain columns between -.5 & 0,5
ii = rho>-.5 & rho < .5;
fprintf('Features Retained\n');
fprintf('%s\n',newHeader{1,ii});
allD = allD(:,ii);
newHeader = newHeader(1,ii);

% Remove Lower 1% and upper 1% data for each control
uControls = unique(allTxt);
jj = false(size(allTxt,1),1);
for i = 1:numel(uControls)
    ii = find(strcmpi(uControls{i,:},allTxt));
    kk = allInten(ii,1)>quantile(allInten(ii,1),.05) & allInten(ii,1)<quantile(allInten(ii,:),.95);
    jj(ii(kk)) = true;
end
allD = allD(jj,:);
allTxt = allTxt(jj,:);
allInten = allInten(jj,:);
allMorRatio = allMorRatio(jj,:);
allMorIntensity = allMorIntensity(jj,:);
allTxtOrg = allTxtOrg(jj,:);
fprintf('# Cells removed by lower & upper quartile %i\n',sum(~jj));
% param.meaninc = mean(allD);
% param.varinc = var(allD);
clear D focus cnt tok iFiles mxRw 


% Feature Selection/Reduction
% newHeader = param.datahdr(1,param.datafeat);
featureReduction = true;
numFeatures = 30;
if(numel(newHeader)<=numFeatures)
    featureReduction = false;
end
if(featureReduction)    
    redFeatures = unsupervisedGreedyFS(allD,numFeatures);
else
    redFeatures = true(1,size(allD,2));
end
fprintf('Features Chosen\n')
fprintf('%s\n',newHeader{1,redFeatures});
% [cf]= princomp(allD(:,redFeatures));
allD = allD(:,redFeatures);
% fprintf('System Paused\n');
% pause;
% Remove
for i = 1:numel(uControls)
    ii = (strcmpi(uControls{i,:},allTxt));
    fprintf('%s\t: %d\n',uControls{i,:},sum(ii));     
end

% Perform clustering high number of centroids -
% Uneven number of samples
k = 5;
distanceTypeForSampling = 'Euclidean';
graphType = 'jaccard';
gps = getGroupIndices(allTxt,unique(allTxt));
uG = unique(gps);

grpCentroids = nan(numel(uG),size(allD,2));
% Get randomly sampled data
vk = false(size(allD,1),1);
mx = 5000;
for i = 1:numel(uG)
    ii = find(gps ==uG(i));
    vk(ii(randperm(numel(ii),mx)),1) = true;
end
[ indx ] = phenograph(allD(vk,:),k,'graphtype',graphType);
uInd = unique(indx);
mCent = nan(numel(uInd),size(allD,2));
numPointsPerCluster = nan(numel(uInd),1);
% mGrp = nan(500,1);
% clsGrp = nan(numel(gps),1);
clusterInten = nan(numel(uInd),1);


for j = 1:numel(uInd)
    kk =  indx == uInd(j);
    mCent(j,:)= mean(allD(kk,:));
    clusterInten(j,1)= mean(allInten(indx == uInd(j),:));
    numPointsPerCluster(j,1) = sum(indx == uInd(j));
end

clear options grpData i j cnt ii kk ;
fprintf('Computed Centroids\n');
clusterInten = (clusterInten - min(clusterInten))./...
    (max(clusterInten) - min(clusterInten));
% Number of clusters Per Group
%% Assign all the data to the nearest centroid
nallLabel = knnclassify(allD, mCent, 1:size(mCent,1));
radiusMarker = zeros(numel(uInd),1);
for i = 1:numel(uInd)
    ii = nallLabel == uInd(i);
    radiusMarker(i) = sum(ii);
    fprintf('%d\n',sum(ii));
end
radiusMarker = (radiusMarker - min(radiusMarker))./...
                    (max(radiusMarker)-min(radiusMarker));
radiusMarker = radiusMarker*100+40;
% Visualize data
pcaRedData = compute_mapping(mCent,'t-SNE',numDimsVis);
%% Number of clusters Per Group

maps = jet(numel(uInd));
figure; hold on;
for i = 1:1
%     ii = mGrp == i;
    if(numDimsVis == 3)
        scatter3(pcaRedData(:,1),pcaRedData(:,2),pcaRedData(:,3),...
            radiusMarker(ii),maps(i,:),'filled');
        
    else
        scatter(pcaRedData(:,1),pcaRedData(:,2),radiusMarker(:),...
            maps(i,:),'filled');
    end
end
hold off;
if(numDimsVis == 3)
    view([az el]);
end
grid on;
% set(gca,'xtick',[]);set(gca,'ytick',[]);set(gca,'ztick',[]);
title('Level-1 Centroids');
legend(uControls);



%% Up sample data - Compute enrichment per cluster - Technique 1
close all;

distMat = zeros(numel(uI2),max(gps));
for i = 1:numel((uI2))
    ii = nallLabel==i;
    for j = 1:max(gps)
        jj = gps == j;
        distMat(i,j) = sum(and(ii,jj));
    end
end

pvalDist = bsxfun(@rdivide,distMat,sum(distMat,1));
pvalDist = pdist2(pvalDist,pvalDist,'euclidean');
% pvalDist(isnan(pvalDist)) = inf;
% pvalDist = pdist2(mCent,mCent,'euclidean');
[~,st,~] = kruskal(1-eye(size(pvalDist,1)),pvalDist);
% Plot lines on the 3D plot
nWeight = zeros(size(st,1),1);
for i = 1:size(st)
    nWeight(i,1) = pvalDist(st(i,1),st(i,2));
end
% Visualize the graph
figure;hold on;
if(numDimsVis == 3)
    scatter3(pcaRedData(:,1),pcaRedData(:,2),pcaRedData(:,3),...
        radiusMarker,maps,'filled');
    for i = 1:size(st,1)
        xv = [pcaRedData(st(i,1),1) pcaRedData(st(i,2),1) ];
        yv = [pcaRedData(st(i,1),2) pcaRedData(st(i,2),2) ];
        zv = [pcaRedData(st(i,1),3) pcaRedData(st(i,2),3) ];
        line(xv,yv,zv,'Color',[.5 .5 .5]);
    end
    view([az el]);
else
    scatter(pcaRedData(:,1),pcaRedData(:,2),radiusMarker,maps,'filled');
    for i = 1:size(st,1)
        xv = [pcaRedData(st(i,1),1) pcaRedData(st(i,2),1) ];
        yv = [pcaRedData(st(i,1),2) pcaRedData(st(i,2),2) ];        
        line(xv,yv,'Color',[.5 .5 .5]);
    end
end
% for i = 1:numel(uI2)
%     if(numDimsVis == 3)
%         text(wPcaRedData(i,1),wPcaRedData(i,2),wPcaRedData(i,3),num2str(i));
%     else
%         text(wPcaRedData(i,1),wPcaRedData(i,2),num2str(i));
%     end
% end
hold off;
grid on;
% xlim(x1);ylim(y1);

% Plot individual landmarks on a subplot

% #Subplots = 15
sy=4;sx = 3;
dm = bsxfun(@rdivide,distMat,sum(distMat));


for i = 1:max(gps)
    figure;
    if(numDimsVis == 3)
        scatter3(pcaRedData(:,1),pcaRedData(:,2),pcaRedData(:,3),...
            radiusMarker,dm(:,i),'filled');
        view([az el]);
    else
        scatter(pcaRedData(:,1),pcaRedData(:,2),...
            radiusMarker,dm(:,i),'filled');
    end
    caxis([0 .5]);
%     set(gca,'Xtick',[]);set(gca,'ytick',[]);set(gca,'ztick',[]);
%     colorbar;
    title(uControls{i,:});
%     axis image;
end
%% END*******************************


%