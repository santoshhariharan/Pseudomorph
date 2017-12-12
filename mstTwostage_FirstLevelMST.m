%% Clear and close everything
clear; clc;close all;
%% Load centroid files and recluster using Phenograph

% Load the centroid file
centroidFileName = 'centroidPerControl_Feat160F_5K.mat';
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
%% Define Color Map
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
hold off; axis off;
% Cluster Map
clsColor = colorcube(numel(uIndx));
figure; hold on;
for i = 1:numel(uIndx)
    ii = clusterIndex == uIndx(i);
    scatter(redData(ii,1),redData(ii,2),mFraction2(ii,1),...
        clsColor(i,:),'filled'); 
end
hold off;axis off;
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
        disMat(i,j) = sum((kk));
    end
end
radiusCls = radiusCls./sum(radiusCls);
pval = enrich_score(disMat');
% radius = (radius - min(radius))./(max(radius)-min(radius));
%% Compute Ensemble MST
mstCutoffValue= 0.5;
numK = floor(.01*size(mCent,1));
[row,col,weight,aPicked] = getGraphByEMST(mCent,numK,mstCutoffValue,100);
fprintf('\n');
%% Plot MST Data
mstCutoffValue= 0.5;
ii = weight>=mstCutoffValue;
r = row(ii);
c = col(ii);
v = weight(ii);
% Get Edges
edgeWidth = (v - min(v))./(max(v)-min(v));
edgeWidth = 1*edgeWidth+.05;
edgeWidthMatrix = zeros(size(mCent,1));
for i = 1:numel(r)
    edgeWidthMatrix(r(i),c(i)) = v(i);
    edgeWidthMatrix(c(i),r(i)) = v(i);
end
% viewMSTPie2(redData(uIndx,:),disMat,maps,[r c],edgeWidth,.05,radiusCls);

figure; hold on;
for i = 1:numel(r)
    coord = redData([r(i);c(i)],:)';
    line(coord(1,:),coord(2,:),'Linewidth',edgeWidth(i),'Color',1-(v(i)*[1 1 1]));
end

for i = 1:max(mGrp)
    ii = mGrp == i;
    scatter(redData(ii,1),redData(ii,2),mFraction2(ii,1),...
        maps(i,:),'filled'); 
end
hold off; axis off;
%% Find shortest path between clusters
% Create Adjacency Matrix
adjacencyMatrix = zeros(size(mCent,1));
distanceMatrix = zeros(size(mCent,1));
for i = 1:numel(r)
    adjacencyMatrix(r(i),c(i)) = 1;
    adjacencyMatrix(c(i),r(i)) = 1;
    distanceMatrix(r(i),c(i)) = pdist2(mCent(r(i),:),mCent(c(i),:));
    distanceMatrix(c(i),r(i)) =distanceMatrix(r(i),c(i)); 
end
cc = cell(200,1);cnt = 1;
rowColumnClusterIndex = zeros(200,1);
% uIndx = sort(uIndx);
for i = 1:numel(uIndx)    
    for j = i+1:numel(uIndx)
        [~, path] = dijkstra(adjacencyMatrix,distanceMatrix,uIndx(i),...
        uIndx(j));
        cc{cnt,1} = path;
        rowColumnClusterIndex(cnt,1) = uIndx(i);
        rowColumnClusterIndex(cnt,2) = uIndx(j);
        cnt=cnt+1;
    end
end
cc = cc(1:cnt-1,1);
rowColumnClusterIndex =rowColumnClusterIndex(1:cnt-1,:);
disp('DONE')
%% Create new MST adjacency matrix
adjacencyMatrixShort = zeros(size(mCent,1));
for i = 1:numel(cc)
    pIndex = cc{i,1};
    for j = 2:numel(pIndex)
        adjacencyMatrixShort(pIndex(j-1),pIndex(j)) = 1;
        adjacencyMatrixShort(pIndex(j),pIndex(j-1)) = 1;        
    end
end
%% Plot Shortest MST data using function
viewMSTPie2(redData,clusterIndex,mGrp,...
                                    maps,adjacencyMatrixShort,edgeWidthMatrix,...
                                    controlNames)
%% Retain Clusters between Groups & Plot Connections
% grp2Retain = [4 5 6 7 8 9 10 11 13];
grp2Retain = [1 2 3 6 7 12 13];
grp2Retain = [ 4 5 6 ];
% Find clusters where the group has atleast 10% of cells
dMat = bsxfun(@rdivide,disMat,sum(disMat,2));
dMat = dMat>.1;
dMat = dMat(:,grp2Retain);
cls2Retain = sum(dMat,2)>0;
adjacencyMatrixGroup = zeros(size(mCent,1));
% grp2Retain = [4 5];
% Pull all clusters that have statistical enrichment
% pvalEnr = pval<=(.05/(numel(uIndx)*max(mGrp)));
% cls2Retain = sum(pvalEnr(grp2Retain,:))>0;
%  cls2Retain = sum(pvalEnr(grp2Retain,:))>0;
cls2Retain = sort(uIndx(cls2Retain));

% Get connection data based on clusters to retain
tmpPath = cell(10,1);
tmpIndex = nan(10,2);
cnt = 1;
for i = 1:numel(cls2Retain)
    ii = rowColumnClusterIndex(:,1) == cls2Retain(i);
    for j = i+1:numel(cls2Retain)
        jj = rowColumnClusterIndex(:,2) == cls2Retain(j);
        kk = logical(ii.*jj);
        tmpPath(cnt,1) = cc(kk,1);
        pIndex = cc{kk,1};
        for k = 2:numel(pIndex)
            adjacencyMatrixGroup(pIndex(k-1),pIndex(k)) = 1;
            adjacencyMatrixGroup(pIndex(k),pIndex(k-1)) = 1;
        end
        tmpIndex(cnt,1) = cls2Retain(i);
        tmpIndex(cnt,2) = cls2Retain(j);
        cnt = cnt+1;
    end
end
tmpPath = tmpPath(1:cnt-1,1);
tmpIndex = tmpIndex(1:cnt-1,:);
%
nClsIndex = false(numel(clusterIndex),1);
for i = 1:numel(cls2Retain)
    ii = clusterIndex==cls2Retain(i);
    nClsIndex(ii) = true;
end
nCindex = clusterIndex;
nCindex(~nClsIndex) = 0;
viewMSTPie2(redData,nCindex,mGrp,...
                                    maps,adjacencyMatrixGroup,edgeWidthMatrix,...
                                    controlNames);
%%
figure; hold on;
% Plot data in gray
% scatter(redData(:,1),redData(:,2),mFraction2,[.8 .8 .8],'filled'); 
% Plot Lines  
for i = 1:numel(tmpPath)
    pIndex = tmpPath{i,1};
    if(isnan(pIndex))
        continue;
    end
    for j = 2:numel(pIndex)
        coord = redData([pIndex(j-1);pIndex(j)],:)';
        line(coord(1,:),coord(2,:),'Linewidth',.5,'Color',[.7 .7 .7]);
    end
end

% Plot the points
for i = 1:numel(tmpPath)
    pIndex = tmpPath{i,1};
    scatter(redData(pIndex,1),redData(pIndex,2),mFraction2(pIndex,1),...
        maps(mGrp(pIndex),:),'filled');
end

% Plot Clusters centers

% Plot Centers
for i = 1:numel(cls2Retain)
    scatter(redData(cls2Retain(i),1),redData(cls2Retain(i),2),100,...
        maps(mGrp(cls2Retain(i)),:),'filled'); 
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
vK = 1;
[idx] = knnsearch(newData,allD,'K',vK);
%% Compute Histogram
uTxt = unique(allTxt);
CP = zeros(numel(uTxt),size(newData,1));
for i = 1:numel(uTxt)
    ii = strcmpi(allTxt,uTxt{i,:});
    for j = 1:size(newData,1)
        jj = idx == j;
        CP(i,j) = sum(ii.*jj);
    end
%     CP(i,:) = histcounts(cellScore(ii),[-.05:.05:1]);    
end
CP = bsxfun(@rdivide,CP,sum(CP,2));
%%
[r2, c2] = find(nAdjMatrix > 0);
figure; hold on;
for i = 1:size(r2,1)  
    coord = redDataTmp([r2(i);c2(i)],:)';
    line(coord(1,:),coord(2,:),'Linewidth',.5,'Color',[.7 .7 .7]);    
end


% % Plot Centers
% for i = 1:numel(uIndx)
%     scatter(redData(uIndx(i),1),redData(uIndx(i),2),100,...
%         maps(mGrp(uIndx(i)),:),'filled'); 
% end


% Plot Histogram
scatter(redDataTmp(:,1),redDataTmp(:,2),mFracTmp,...
        CP(3,:)','filled');
    hold off;
return;
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
