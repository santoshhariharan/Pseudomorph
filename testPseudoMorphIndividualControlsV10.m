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
pth='F:\Projects\Proteinlocalization\PseudoMorph\DataFiles';
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
% clustersPerLandmark = true;

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
% Remove nan entries
ii = sum(isnan(allD),2) ==0;
allD = allD(ii,:);
allInten = allInten(ii,:);
allTxt = allTxt(ii,:);
allMorRatio = allMorRatio(ii,:);
allMorIntensity = allMorIntensity(ii,:);
allTxtOrg= allTxtOrg(ii,:);
fprintf('Removed NAN\n');

% Remove incorrectly segmented & low intensity Objects
ii = allInten>30 & allMorRatio <= .7;
allD = allD(ii,:);
allTxt = allTxt(ii,:);
allTxtOrg = allTxtOrg(ii,:);
allInten = allInten(ii,:);
allMorRatio = allMorRatio(ii,:);
allMorIntensity = allMorIntensity(ii,:);
fprintf('Removed %i of %i\n',sum(~ii),numel(ii));

% Normalization type
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
param.meaninc = mean(allD);
param.varinc = var(allD);


% Feature Selection/Reduction
% newHeader = param.datahdr(1,param.datafeat);
featureReduction = true;
numFeatures = 10;
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
% Print number of cells per control
for i = 1:numel(uControls)
    ii = (strcmpi(uControls{i,:},allTxt));
    fprintf('%s\t: %d\n',uControls{i,:},sum(ii));    
end


clear D focus cnt tok iFiles mxRw 
clear allMorRatio allInten
clear ii jj kk rho mxRw 
clear intensityFeature intFeat roiIntFeat nucAreaFeat cellAreaFeat
clear fNames filePrefix randPrc
% Perform clustering high number of centroids -
% Uneven number of samples
k = 5;
distanceTypeForSampling = 'Euclidean';
graphType = 'jaccard';
gps = getGroupIndices(allTxt,unique(allTxt));
uG = unique(gps);
mCent = nan(500,size(allD,2));
mGrp = nan(500,1);
clsGrp = nan(numel(gps),1);
clusterInten = nan(500,1);
numPointsPerCluster = nan(500,1);
grpCentroids = nan(numel(uG),size(allD,2));
cnt = 1;
for i = 1:numel(uG)
    ii = find(gps ==uG(i));
    grpData = allD(ii,:);
    grpInten = allMorIntensity(ii,1);
    grpCentroids(i,:) = mean(grpData);
%     jj = getDensityBasedSampling(grpData,distanceTypeForSampling
    jj  = true(numel(ii),1);
    [ indx ] = phenograph(grpData(jj,:),k,'graphtype',graphType);
    uInd = unique(indx);
    for j = 1:numel(uInd)
       kk =  indx == uInd(j);
       clsGrp(ii(kk)) = cnt;
       mCent(cnt,:)= mean(grpData(kk,:));
       clusterInten(cnt,1)= mean(grpInten(indx == uInd(j),:));
       mGrp(cnt,1) = i;
       numPointsPerCluster(cnt,1) = sum(indx == uInd(j));
       cnt = cnt+1;
    end     
end
%
ii = sum(isnan(mCent),2)==0;
mCent = mCent(ii,:);
mGrp = mGrp(ii,1);
clusterInten = clusterInten(ii,1);
numPointsPerCluster = numPointsPerCluster(ii,1);
% Compute Level 2 centroids
[indxN,~] = phenograph(mCent,k);
fprintf('Computed Centroids\n');
clusterInten = (clusterInten - min(clusterInten))./...
    (max(clusterInten) - min(clusterInten));

clear options grpData i j uInd indx cnt ii;
clear comm kk grpInten
% Number of clusters Per Group
radiusMarker = zeros(numel(mGrp),1);
for i = 1:numel(uControls)
    radiusMarker(mGrp==i) = numPointsPerCluster(mGrp==i)./...
                            sum(numPointsPerCluster(mGrp==i));
    fprintf('%d\n',sum(mGrp==i));
end
radiusMarker = (radiusMarker - min(radiusMarker))./...
                    (max(radiusMarker) - min(radiusMarker));
radiusMarker = radiusMarker*80 + 20;
% Visualize data
pcaRedData = compute_mapping(mCent,'t-SNE',numDimsVis);
%% Reassign centroids to Group
% Compute prior probability
% Minimum number of centroids per control
minPoints = 500;
for i = 1:numel(uControls)    
    fprintf('%d\n',sum(mGrp==i));
    minPoints = min(minPoints,sum(mGrp==i));
end
kk = floor(minPoints/2);
priorProb = zeros(max(mGrp),1);
for i = 1:max(mGrp)
    priorProb(i) = sum(mGrp==i);
end
priorProb = priorProb-1;
priorProb = priorProb./sum(priorProb);
idx = knnsearch(mCent,mCent,'k',kk+1);
idx = idx(:,2:end);
idx = mGrp(idx);
idxn = hist(idx',[1:max(gps)])';
[~,maxIdxn] = max(idxn,[],2);
xx = mGrp ~= maxIdxn;
idxn = idxn(xx,:);
pvalue = nan(size(idxn,1),1);
pvalue2 = nan(size(idxn,1),1);
pp = [find(xx) mGrp(xx) maxIdxn(xx)];
for i = 1:size(idxn,1)
    m= idxn(i,pp(i,3));
    m2 = idxn(i,pp(i,2));
    if(m == 1 || m2 == 1)
        continue;
    end    
    pvalue(i,1) = 1-binocdf(m-1,kk,priorProb(pp(i,3)));
    pvalue2(i,1) = 1-binocdf(m2-1,kk,priorProb(pp(i,2)));
end
ii = isnan(pvalue) | isnan(pvalue2);
pvalue = pvalue(~ii);
pvalue2 = pvalue2(~ii);
xx = (pvalue < (.05./(numel(mGrp).*max(mGrp))));
yy = (pvalue2 >= (.05./(numel(mGrp).*max(mGrp))));
xx = and(xx,yy);
pp = pp(xx,:);

% Plot matrix
mergeMat = zeros(max(mGrp),max(mGrp));
for i = 1:max(mGrp)
    ii = pp(:,2) == i;
    for j = 1:max(mGrp)
        jj = pp(:,3) == j;
        mergeMat(i,j) = sum(ii.*jj);
    end
end

newMGrp = gps;
for i = 1:size(pp,1)
    ii = clsGrp == pp(i,1);
    newMGrp(ii) = pp(i,3);
end
% newMGrp(xx,1) = pp(:,3);
% pvalue = pvalue < (.05/max(mGrp));
% pp = [mGrp(xx) maxIdxn(xx) pvalue(:,maxIdxn(xx))];

clear minPoints kk priorProb xx yy idx idxn
clear pvalue pvalue2 maxidxn
%% Number of clusters Per Group

figure;
if(numDimsVis == 3)
    scatter3(pcaRedData(:,1),pcaRedData(:,2),pcaRedData(:,3),...
        radiusMarker,clusterInten,'filled');
    [az, el] = view;
else
    scatter(pcaRedData(:,1),pcaRedData(:,2),...
        radiusMarker,clusterInten,'filled');    
end
set(gca,'xtick',[]);set(gca,'ytick',[]);set(gca,'ztick',[]);
figure; hold on;
for i = 1:numel(uG)
    ii = mGrp == i;
    if(numDimsVis == 3)
        scatter3(pcaRedData(ii,1),pcaRedData(ii,2),pcaRedData(ii,3),...
            radiusMarker(ii),maps(i,:),'filled');
        
    else
        scatter(pcaRedData(ii,1),pcaRedData(ii,2),radiusMarker(ii),...
            maps(i,:),'filled');
    end
end
hold off;
if(numDimsVis == 3)
    view([az el]);
end
grid on;
set(gca,'xtick',[]);set(gca,'ytick',[]);set(gca,'ztick',[]);
title('Level-1 Centroids');
legend(uControls);

uI2 = unique(indxN);
nMP = jet(numel(uI2));
figure; hold on;
for i = 1:numel(uI2)
    ii = indxN == uI2(i);
    if(numDimsVis == 3)
        scatter3(pcaRedData(ii,1),pcaRedData(ii,2),pcaRedData(ii,3),...
            radiusMarker(ii),nMP(i,:),'filled');
    else
        scatter(pcaRedData(ii,1),pcaRedData(ii,2),radiusMarker(ii),...
            nMP(i,:),'filled');
    end
end
hold off;
if(numDimsVis == 3)
    view([az el]);
end
title('Level-2 Centroids');
set(gca,'xtick',[]);set(gca,'ytick',[]);set(gca,'ztick',[]);
grid on;
x1 = xlim;
y1 = ylim;
 
%% Compute location of nearest centroids

wPcaRedData = zeros(numel(uI2),numDimsVis);
wPcaRedDataAllD = zeros(numel(uI2),size(allD,2));
distData = zeros(numel(uI2),max(mGrp));
nRadMarker = zeros(numel(uI2),1);
for i = 1:numel(uI2)
    ii = indxN == uI2(i);
    nMGrp = mGrp(ii);
    distData(i,:) = hist(mGrp(ii),[1:max(mGrp)]);
    nPcaRedData=bsxfun(@times,pcaRedData(ii,:),numPointsPerCluster(ii));
    nPcaRedData = nPcaRedData./sum(numPointsPerCluster(ii));
    wPcaRedData(i,:) = sum(nPcaRedData);
    nRadMarker(i) = sum(numPointsPerCluster(ii));
    
    nPcaRedData=bsxfun(@times,mCent(ii,:),numPointsPerCluster(ii));
    nPcaRedData = nPcaRedData./sum(numPointsPerCluster(ii));
    wPcaRedDataAllD(i,:) = sum(nPcaRedData);
end
nRadMarker = (nRadMarker./sum(nRadMarker));
nRadMarker = (nRadMarker - min(nRadMarker))./(max(nRadMarker)-min(nRadMarker));
nRadMarker = 160*nRadMarker+40;
figure;
if(numDimsVis == 3)    
    scatter3(wPcaRedData(:,1),wPcaRedData(:,2),wPcaRedData(:,3),...
        nRadMarker,nMP,'filled');
    view([az el]);
else
    scatter(wPcaRedData(:,1),wPcaRedData(:,2),nRadMarker,nMP,'filled');
end
title('New Centroids');
grid on;
xlim(x1);ylim(y1);



%% Up sample data - Compute enrichment per cluster - Technique 1
close all;
% Assign all the data to the nearest centroid
nallLabel = knnclassify(allD, wPcaRedDataAllD, 1:numel(uI2));
distMat = zeros(numel(uI2),max(gps));
for i = 1:numel((uI2))
    ii = nallLabel==i;
    for j = 1:max(gps)
        jj = newMGrp == j;
        distMat(i,j) = sum(and(ii,jj));
    end
end
% ii = sum(distMat(:,2),2)>0;
% distMat = distMat(ii,:);
% wPcaRedDataAllD = wPcaRedDataAllD(ii,:);
% wPcaRedData = wPcaRedData(ii,:);
% pvalDist = enrich_score(distMat);
% pvalDist = (pvalDist< (.05./numel(pvalDist)));

pvalDist = bsxfun(@rdivide,distMat,sum(distMat,1));
newX = compute_mapping(pvalDist,'t-SNE',2);
pvalDist = pdist2(pvalDist,pvalDist,'correlation');
% pvalDist(isnan(pvalDist)) = inf;
% pvalDist = pdist2(wPcaRedDataAllD,wPcaRedDataAllD,'euclidean');
[~,st,~] = kruskal(1-eye(size(pvalDist,1)),pvalDist);
% Plot lines on the 3D plot
nWeight = zeros(size(st,1),1);
for i = 1:size(st)
    nWeight(i,1) = pvalDist(st(i,1),st(i,2));
end
% Visualize the graph
figure;hold on;
if(numDimsVis == 3)
    scatter3(wPcaRedData(:,1),wPcaRedData(:,2),wPcaRedData(:,3),...
        nRadMarker,nMP,'filled');
    for i = 1:size(st,1)
        xv = [wPcaRedData(st(i,1),1) wPcaRedData(st(i,2),1) ];
        yv = [wPcaRedData(st(i,1),2) wPcaRedData(st(i,2),2) ];
        zv = [wPcaRedData(st(i,1),3) wPcaRedData(st(i,2),3) ];
        line(xv,yv,zv,'Color',[.5 .5 .5]);
    end
    view([az el]);
else
%     scatter(wPcaRedData(:,1),wPcaRedData(:,2),nRadMarker,nMP,'filled');
    scatter(wPcaRedData(:,1),wPcaRedData(:,2),nRadMarker,nMP,'filled');
    for i = 1:size(st,1)
        xv = [wPcaRedData(st(i,1),1) wPcaRedData(st(i,2),1) ];
        yv = [wPcaRedData(st(i,1),2) wPcaRedData(st(i,2),2) ];        
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
%%
viewScatterPie2(wPcaRedData,nallLabel,allTxt,maps);
% xlim(x1);ylim(y1);

% Plot individual landmarks on a subplot

% #Subplots = 15
% sy=4;sx = 3;
% dm = bsxfun(@rdivide,distMat,sum(distMat));
% 
% 
% for i = 1:size(maps,1)
%     figure;
%     if(numDimsVis == 3)
%         scatter3(wPcaRedData(:,1),wPcaRedData(:,2),wPcaRedData(:,3),...
%             nRadMarker,dm(:,i),'filled');
%         view([az el]);
%     else
%         scatter(wPcaRedData(:,1),wPcaRedData(:,2),...
%             nRadMarker,dm(:,i),'filled');
%     end
%     caxis([0 .5]);
% %     set(gca,'Xtick',[]);set(gca,'ytick',[]);set(gca,'ztick',[]);
% %     colorbar;
%     title(uControls{i,:});
% %     axis image;
% end
%% END*******************************


%% Find distances between centroids
% Find candidates for cluster merging by controls
[closestDisIdx,closestDis] = knnsearch(mCent,mCent,'k',2);
closestDisIdx = closestDisIdx(:,2:end);

% Remove points whose closest two neighbours are from same group
a = mGrp(closestDisIdx(:,1));
% b = mGrp(closestDisIdx(:,2));
iiIndex = find(a ~= mGrp); % Find the ones which are closest to other groups
closestDisIdx = closestDisIdx(iiIndex,:);
closestDis = closestDis(iiIndex,2);

% Remove duplicates --
xR = false(numel(iiIndex),1);
uCIndex = unique(closestDisIdx);
for i = 1:numel(uCIndex)
    ii = find(closestDisIdx == uCIndex(i));    
    [~,p] = min(closestDis(ii));
    xR(ii(p)) = true;
end
closestDisIdx = closestDisIdx(xR);
closestDis = closestDis(xR);
iiIndex = iiIndex(xR);

% Plot potential candidates
figure; hold on;
for i = 1:numel(uG)
    ii = mGrp == i;
    scatter3(pcaRedData(ii,1),pcaRedData(ii,2),pcaRedData(ii,3),radiusMarker(ii),...
                    [.5 .5 .5],'filled');
end
nMap = jet(numel(iiIndex));
for i = 1:numel(iiIndex)
    scatter3(pcaRedData(iiIndex(i),1),pcaRedData(iiIndex(i),2),pcaRedData(iiIndex(i),3),...
        radiusMarker(iiIndex(i)),...
                    nMap(i,:),'filled');
    scatter3(pcaRedData(closestDisIdx(i,1),1),pcaRedData(closestDisIdx(i,1),2),...
        pcaRedData(closestDisIdx(i,1),3),...
        radiusMarker(closestDisIdx(i,1)),...
                    nMap(i,:),'filled');            
end
hold off;
view([az el]);grid on;
title('Centroids for Potential Merging');
%% Assign cluster labels for all data based on each group
allLbl = zeros(numel(gps),1);
uGrp = unique(gps);
for i = 1:numel(uGrp)
    ii = find(mGrp == uGrp(i));
    jj = gps == uGrp(i);
    allLbl(jj,1)  = knnclassify(allD(jj,:), mCent(ii,:), ii);
end

% Check if clusters need to be merged

allAccu = zeros(numel(iiIndex),3);
nPerm = 100;
fprintf('Checking Mergers..............');
for i = 1:numel(iiIndex)      
    fprintf('\b\b\b\b\b\b\b\b%7.3f%%',i*100/numel(iiIndex));
    [allAccu(i,1), allAccu(i,2), allAccu(i,3)] = ...
                        permByClassification( allD(allLbl==iiIndex(i),:),...
                            allD(allLbl==closestDisIdx(i),:),nPerm );
end
fprintf('\n');

%% Plot errors on the same t-SNE plot
allErr = 1- allAccu(:,2);
figure; hold on;
for i = 1:numel(uG)
    ii = mGrp == i;
    scatter(pcaRedData(ii,1),pcaRedData(ii,2),radiusMarker(ii),...
                    [.5 .5 .5],'filled');
end
scatter(pcaRedData(iiIndex,1),pcaRedData(iiIndex,2),...
                        radiusMarker(iiIndex),allErr,'filled');
scatter(pcaRedData(closestDisIdx,1),pcaRedData(closestDisIdx,2),...
                        radiusMarker(closestDisIdx),allErr,'filled');       
hold off;                    
%% Merge clusters where error rate is greater than .2
% Index of clusters to be merged
nLbl = allLbl;
ii = find(allErr>.2);

% Assuming merging driven by larger centroid
for i = 1:numel(ii)
    jj = nLbl == iiIndex(ii(i));
    kk = nLbl == closestDisIdx(ii(i));
    if(sum(jj)>=sum(kk))
        nLbl(kk) = iiIndex(ii(i));
    else
        nLbl(jj) = closestDisIdx(ii(i));
    end    
end

nULbl = unique(nLbl);
nMCent = zeros(numel(nULbl),size(allD,2));
% nMGrp = zeros(numel(nULbl),1);
for i= 1:numel(nULbl)
    nMCent(i,:) = mean(allD(nLbl == nULbl(i),:));
%     nMGrp(i,1) = mGrp(nULbl(i));
end

% Assign all the data -> Upsample

nLbl = knnclassify(allD, nMCent, 1:numel(nULbl));
nPcaRedData = compute_mapping(nMCent,'t-SNE',2);
viewScatterPie2(nPcaRedData,nLbl,allTxt,maps);
viewScatterPie2(nPcaRedData,gps,allTxt,maps);
