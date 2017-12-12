
%
clear; clc;
% Load the centroid file
centroidFileName = 'centroidPerControl_Feat49F_5K.mat';
unKnowmRootDir = 'F:\Projects\Proteinlocalization\PseudoMorph\Unknowns';

load(centroidFileName);
% Cluster and compute centroids
reclusterInitialCentroids = true;

if(reclusterInitialCentroids)
    k = 5;
    [clusterIndex,mCent2,mFrac2,uIndx] = clusterByPhenograph( mCent,k );
end
%% Compute dimension reduction
redData = compute_mapping(mCent,'t-SNE',2);
%% Define Map
load('maps.mat');
minSize = 20;
maxSize = 70;
mFraction2 = maxSize*mFraction+minSize;
%% Plot data - Centers on low dimensional plot
figure; hold on;
for i = 1:max(mGrp)
    ii = mGrp == i;
    scatter(redData(ii,1),redData(ii,2),mFraction2(ii,1),...
        maps(i,:),'filled'); 
end
hold off;  
% Cluster Map
clsColor = colorcube(numel(uIndx));
figure; hold on;
for i = 1:numel(uIndx)
    ii = clusterIndex == uIndx(i);
    scatter(redData(ii,1),redData(ii,2),mFraction2(ii,1),...
        clsColor(i,:),'filled'); 
end
hold off;
%% Compute cluster Distribution matrix
% uIndx = unique(clusterIndex);
disMat = zeros(numel(uIndx),max(mGrp));
radiusCls = zeros(numel(uIndx),1);
for i = 1:numel(uIndx)
    ii = clusterIndex==uIndx(i);
    radiusCls(i) = sum(ii);
    for j = 1:max(mGrp)
        jj = mGrp==j;
        kk = logical(ii.*jj);
        disMat(i,j) = sum(mFraction(kk));
    end
end
radiusCls = radiusCls./sum(radiusCls);
% radius = (radius - min(radius))./(max(radius)-min(radius));
%% Compute Ensemble MST
mstCutoffValue= .5;
numK = size(mCent2,1);
[r,c,v] = getGraphByEMST( mCent2,numK,mstCutoffValue,500 );
%% View data by Pie chart
edgeWidth = (v - min(v))./(max(v)-min(v));
edgeWidth = 1*edgeWidth+.1;
viewMSTPie2(redData(uIndx,:),disMat,maps,[r c],edgeWidth,.05,radiusCls);
%% Create a KNN graph per pair of clusters
adjacencyMatrix = zeros(size(mCent,1));
numClusterObjects = mFrac2*numel(clusterIndex);
numK = ceil(quantile(numClusterObjects,.1));
for i = 1:numel(r)
    ii = find(clusterIndex == uIndx(r(i)));
    jj = find(clusterIndex == uIndx(c(i)));  
    allInd = [ii;jj];
    clustePairData = [mCent(ii,:);mCent(jj,:)];
%     clusterCenterIndexInClusterPair = [find(allInd==uIndx(r(i)));find(allInd==uIndx(c(i)))];
    m = size(clustePairData,1);
    [idx,dis] = knnsearch(clustePairData,clustePairData,'K',numK+1);
    idx = idx(:,2:end);
    dis = dis(:,2:end);
    for j = 1:m
        adjacencyMatrix(allInd(j),allInd(idx(j,:))) = 1;
    end
    
    adjacencyMatrix = adjacencyMatrix+adjacencyMatrix';
    adjacencyMatrix = adjacencyMatrix>0;    
%     allInd(idx)
end
% Remove Straight Connections
adjacencyMatrix(uIndx,uIndx) = false;
[r1,c1] = find(adjacencyMatrix);
%% Plot KNN Graph per pair of clusters
% Plot data

figure; hold on;
for i = 1:size(r1,1)  
    coord = redData([r1(i);c1(i)],:)';
    line(coord(1,:),coord(2,:),'Linewidth',.5,'Color',[.7 .7 .7]);    
end

for i = 1:max(mGrp)
    ii = mGrp == i;
    scatter(redData(ii,1),redData(ii,2),mFraction2(ii,1),...
        maps(i,:),'filled'); 
end
hold off;  
%% Compute shortest Path Across clusters
distanceMatrix = pdist2(mCent,mCent);
distanceMatrix(~adjacencyMatrix) = inf;
cc = cell(200,1);cnt = 1;
for i = 1:numel(r)
    [~, path] = dijkstra(adjacencyMatrix,distanceMatrix,uIndx(r(i)),...
        uIndx(c(i)));
    cc{cnt,1} = path;
    cnt=cnt+1;
end
cc = cc(1:cnt-1,1);
%% Plot MST data
figure; hold on;
% Plot data in gray
% scatter(redData(:,1),redData(:,2),mFraction2,[.8 .8 .8],'filled'); 
% Plot Lines  
tmp = [];
for i = 1:numel(cc)
    pIndex = cc{i,1};
    if(isnan(pIndex))
        continue;
    end
    tmp = [tmp;pIndex'];
    for j = 2:numel(pIndex)
        coord = redData([pIndex(j-1);pIndex(j)],:)';
        line(coord(1,:),coord(2,:),'Linewidth',.5,'Color',[.7 .7 .7]);
    end
end

% Plot Centers
for i = 1:numel(uIndx)
    scatter(redData(uIndx(i),1),redData(uIndx(i),2),80,...
        maps(mGrp(uIndx(i)),:),'filled'); 
end
tmp = unique(tmp);
mGrpTmp = mGrp(tmp);
uGrp = unique(mGrpTmp);
mFracTmp = mFraction2(mGrpTmp);
redDataTmp = redData(tmp,:);
% Plot Joining points
for i = 1:numel(uGrp)
    ii = mGrpTmp == uGrp(i);
    scatter(redDataTmp(ii,1),redDataTmp(ii,2),mFracTmp(ii,1),...
        maps(uGrp(i),:),'filled'); 
end
hold off;
%% Create New Adjacency Matrix 
% Find the cluster centers in the new matrix
clusCenters = false(numel(tmp,1));
for i = 1:numel(uIndx)
    ii = tmp == uIndx(i);
    clusCenters = or(clusCenters,ii);
end
newData = mCent(tmp,:);
nAdjMatrix = zeros(size(mCent,1));
for i = 1:numel(cc)
    pIndex = cc{i,1};
    if(isnan(pIndex))
        continue;
    end
    for j = 2:numel(pIndex)
        coord = redData([pIndex(j-1);pIndex(j)],:)';
        nAdjMatrix(pIndex(j-1),pIndex(j)) = 1;
    end
end
nAdjMatrix = nAdjMatrix+nAdjMatrix';
nAdjMatrix = nAdjMatrix>0;
nAdjMatrix = nAdjMatrix(tmp,tmp);
nDistanceMatrix = pdist2(newData,newData);
return;
%% Load unknowns with filters

fprintf('Starting pseudomorph\n');
% pth2paramfile='F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data';
load(fullfile(unKnowmRootDir,'parameters.mat'));% Load parameter file
param.rootpath = unKnowmRootDir;
intensityFeature = 'Ch2_INT_Cell_intensity';
intFeat = strcmpi(intensityFeature,param.datahdr);
nucAreaFeat = strcmpi('Ch1_MOR_Nucleus_area',param.datahdr);
cellAreaFeat = strcmpi('Ch1_MOR_Cytoplasm_area',param.datahdr);
filePrefix = '.txt';
fNames = dir(unKnowmRootDir);
columnForControls = 9;
columnForOrganelle = 10;
% featureReduction = true;
% clustersPerLandmark = true;

maxMinTType = false;
% Module 3: Read & Load data after filtering
fprintf('Module 3.......\n');
mxRw = 1000000; 
allD = zeros(mxRw,sum(dataFeat));
allInten = zeros(mxRw,1);
allMorRatio = zeros(mxRw,1);
allTxt = cell(mxRw,1);
allTxtOrg = cell(mxRw,1);
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
    %     Remove cells with spots
    ii = (D.data(:,strcmpi('Ch1_INT_Nucleus_intensity',param.datahdr))./...
        D.data(:,strcmpi('Ch1_INT_Cytoplasm_intensity',param.datahdr)))>3.5;
    jj = (D.data(:,strcmpi('Ch1_INT_Nucleus_intensity_stddev',param.datahdr))./...
        D.data(:,strcmpi('Ch1_INT_Cytoplasm_intensity_stddev',param.datahdr)))>3.5;
    ii = and(ii,jj);      
    D.data = D.data(ii,:);
    D.textdata = D.textdata(ii,:);
    %     Remove cells that have inappropriate nucelar cell ratio
    jj = D.data(:,nucAreaFeat)./(D.data(:,cellAreaFeat)+D.data(:,nucAreaFeat));  
    jj = jj <=.5 & jj >=.2;
    D.data = D.data(jj,:);
    D.textdata = D.textdata(jj,:);
    %     Remove NO intensity cells
    jj = D.data(:,intFeat) >100;
    D.data = D.data(jj,:);
    D.textdata = D.textdata(jj,:);  
    allD(cnt:cnt+size(D.data,1)-1,:) = D.data(:,dataFeat);
    allTxt(cnt:cnt+size(D.data,1)-1,:)= D.textdata(:,columnForControls);
%     allTxt(cnt:cnt+size(D.data,1)-1,:)= D.textdata(:,:);
    allTxtOrg(cnt:cnt+size(D.data,1)-1,:)= D.textdata(:,columnForOrganelle);
    cnt = cnt+size(D.data,1);
end
if(cnt<mxRw)
    allD = allD(1:cnt-1,:);
    allTxt = allTxt(1:cnt-1,:);
    allTxtOrg = allTxtOrg(1:cnt-1,:);
end

fprintf('\n');
% Remove nan entries
ii = sum(isnan(allD),2) ==0;
allD = allD(ii,:);
allTxt = allTxt(ii,:);
allTxtOrg= allTxtOrg(ii,:);
fprintf('\nRemoved NAN %i\n',sum(~ii));

clear D focus cnt tok iFiles mxRw 
clear allMorRatio
clear ii jj kk rho mxRw 
clear intensityFeature intFeat roiIntFeat nucAreaFeat cellAreaFeat
clear fNames filePrefix randPrc

allD = bsxfun(@minus,allD,meanD);
allD = bsxfun(@rdivide,allD,stdD);
%% Compute Nearest Neighbour of Unknowns
vK = 5;
[idx] = knnsearch(newData,allD,'K',vK);
%% Compute Graph distance
[graphDistance] = dijkstra(nAdjMatrix,nDistanceMatrix);
%% Compute cell profiles
cProfile = 0;
for i = 1:vK
    cProfile = cProfile+graphDistance(idx(:,i),clusCenters);
end
cProfile = cProfile/vK;  
cProfile = zscore(cProfile);
cProfileWeights= bsxfun(@rdivide,cProfile,sum(cProfile));
%% Find Nearest Centers By cell Profiles
[~,idx] = min(cProfile,[],2); 
%% Score Each Cell
cellScore = sum(cProfileWeights.*cProfile,2)/size(cProfile,2);
cellScore = (cellScore - min(cellScore))./(max(cellScore) - min(cellScore));
%%
uTxt = unique(allTxt);
CP = zeros(numel(uTxt),numel(0:.05:1));
for i = 1:numel(uTxt)
    ii = strcmpi(allTxt,uTxt{i,:});
    CP(i,:) = histcounts(cellScore(ii),[-.05:.05:1]);    
end
CP = bsxfun(@rdivide,CP,sum(CP,2));
%% Compute Clone Profiles

uTxt = unique(allTxt);
cloneProfile = zeros(numel(uTxt),sum(clusCenters));
for i = 1:numel(uTxt)
    ii = strcmpi(allTxt,uTxt{i,:});
%     cloneProfile(i,:) = histcounts(idx(ii),[.5:1:sum(clusCenters)+.5]);
    cloneProfile(i,:) = sum(cProfile(ii,:));
end
sumCP = sum(cloneProfile,2);
minCP = min(cloneProfile,[],2);
maxCP = max(cloneProfile,[],2);
minMaxCP = maxCP-minCP;
CP = bsxfun(@minus,cloneProfile,minCP);
CP = bsxfun(@rdivide,CP,minMaxCP);
% CP = bsxfun(@rdivide,cloneProfile,sum(cloneProfile,2));
CP = cumsum(CP,2);
% CP = bsxfun(@rdivide,cloneProfile,sumCP);
%% Plot profiles
figure; hold on;
mp = colorcube(numel(uTxt));
for i = 1:numel(uTxt)
    plot(1:size(CP,2),CP(i,:),'-','linewidth',1.5);
end
hold off;
legend(uTxt);
