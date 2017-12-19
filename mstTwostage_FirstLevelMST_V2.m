%% Clear and close everything
clear; clc;close all;
%% Load centroid files 
% Load the centroid file
centroidFileNames = {
%         'F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data\Mito-ER\MITO-MET-PSS1-ER_centroidPerControl_Feat160F_5K.mat'
%         'Sec_centroidPerControl_Feat160F_5K_NoOttawa.mat'
        'F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data\SecPathway\OnlySecPway_centroidPerControl_Feat160F_5K.mat'
%     'F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data\ER\ER_centroidPerControl_Feat160F_5K.mat';
%      'F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data\ERGIC\ERGiC_centroidPerControl_Feat160F_5K.mat';
%      'F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data\Golgi\Golgi_centroidPerControl_Feat160F_5K.mat';               
%      'F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data\Mito-Mam\Mito-MAM_centroidPerControl_Feat160F_5K.mat';               
%      'F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data\Peroxisome-LaminA\Per-LaminA_centroidPerControl_Feat160F_5K.mat';
%      'F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data\Rabs\Rab_centroidPerControl_Feat160F_5K.mat';
%      'F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data\SVPM-Lamp1\SVPM_centroidPerControl_Feat160F_5K.mat' 
     };
unKnowmRootDir = 'F:\Projects\Proteinlocalization\PseudoMorph\Unknowns';
allCent = nan(5000,160);
allTxt = cell(5000,17);
allFrac = nan(5000,1);
cnt = 1;
for i = 1:numel(centroidFileNames)
    load(centroidFileNames{i,:});
    m = size(mCent,1);
    allCent(cnt:cnt+m-1,:) = mCent;
    allTxt(cnt:cnt+m-1,:) = mText;
    allFrac(cnt:cnt+m-1,:) = mFraction;
    cnt = cnt+m;
end
ii = sum(isnan(allCent),2) == 0;
allCent = allCent(ii,:);
allTxt = allTxt(ii,:);
allFrac = allFrac(ii,:);


mCent = allCent;
mFraction = allFrac;
mText = allTxt;

clear allCent allTxt allFrac

controlNames = unique(mText(:,9));
mGrp = getGroupIndices(mText(:,9),controlNames);
load('maps.mat');
% load(centroidFileName);
%% Remove controls
grpNumber = [];
% grpNumber = [8 10 11 12 13 14];
if(~isempty(grpNumber))
    ii = false(numel(mGrp),1);
    mNames = false(numel(controlNames),1);
    for i = 1:numel(grpNumber)
        ii = or(ii,mGrp==grpNumber(i));
        mNames(grpNumber(i)) = true;
    end    
%     mGrp  = mGrp(~ii,:);
    mCent = mCent(~ii,:);
    mFraction = mFraction(~ii,:);
    mText = mText(~ii,:);
    controlNames = (controlNames(~mNames,:));
    maps = maps(~mNames,:);
    mGrp = getGroupIndices(mText(:,9),controlNames);
    % Reset mGrp
%     nMGrp = zeros(numel(mGrp),1);cnt = 1;
%     for i = 1:max(mGrp)
%         ii = mGrp == i;
%         if(sum(ii)>0)
%             nMGrp(ii) = cnt;
%             cnt = cnt+1;
%         end
%     end
%     mGrp = nMGrp;
%     clear nMGrp
end
%% Cluster Using Phenograph
% mCent = zscore(mCent);
clusterCutoffNumPoints = 0;
clusterCutoffNumPoints = clusterCutoffNumPoints*size(mCent,1);
% Cluster and compute centroids
reclusterInitialCentroids = true;

if(reclusterInitialCentroids)
    k = 5;
    [clusterIndex,mCent2,mFrac2,uIndx] = clusterByPhenograph( mCent,k );
end
keepPoints = true(numel(clusterIndex),1);
for i = 1:numel(uIndx)
    ii = clusterIndex == uIndx(i);
    if(sum(ii)<clusterCutoffNumPoints)
        keepPoints(ii) = false;
    end
    
end
mCent = mCent(keepPoints,:);
mGrp  = mGrp(keepPoints,:);
mText = mText(keepPoints,:);
mFraction = mFraction(keepPoints,:);
clusterIndex = clusterIndex(keepPoints,:);
uIndx = unique(clusterIndex);
%% Compute dimension reduction
redData = compute_mapping(mCent,'t-SNE',2);
%% Define Color Map

minSize = 20;
maxSize = 70;
mFraction2 = maxSize*mFraction+minSize;
%% Plot data - Centers on low dimensional plot
figure; hold on;
for i = 1:max(mGrp)
    ii = mGrp == i;
    if(sum(ii)==0)
        continue;
    end
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
%% Compute MST Data
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
%% Plot MST Data
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
costValue = nan(200,1);
rowColumnClusterIndex = zeros(200,3);
% uIndx = sort(uIndx);
for i = 1:numel(uIndx)    
    for j = i+1:numel(uIndx)
        [costValue(cnt), path] = dijkstra(adjacencyMatrix,distanceMatrix,uIndx(i),...
        uIndx(j));
        cc{cnt,1} = path;
        rowColumnClusterIndex(cnt,1) = uIndx(i);
        rowColumnClusterIndex(cnt,2) = uIndx(j);
        rowColumnClusterIndex(cnt,3) = numel(intersect(path,uIndx))-2;
        cnt=cnt+1;
    end
end
cc = cc(1:cnt-1,1);
costValue = costValue(1:cnt-1,1);
rowColumnClusterIndex =rowColumnClusterIndex(1:cnt-1,:);
disp('DONE');
%% Create new MST adjacency matrix
adjacencyMatrixShort = zeros(size(mCent,1));
for i = 1:numel(cc)
    pIndex = cc{i,1};
    for j = 2:numel(pIndex)
        adjacencyMatrixShort(pIndex(j-1),pIndex(j)) = 1;
        adjacencyMatrixShort(pIndex(j),pIndex(j-1)) = 1;        
    end
end
%% Plot Shortest Path between clusters data using function
viewMSTPie2(redData,clusterIndex,mGrp,...
                                    maps,adjacencyMatrixShort,edgeWidthMatrix,...
                                    controlNames)
%% Write Cytoscape file       

writeCytoscapeFile(fileprefix,adjacencyMatrix,edgeWeights,edgeAttribute,...
    clusterIndex,mGrp)
%% Retain Clusters between Groups & Plot Connections
% grp2Retain = [4 5 6 7 8 9 10 11 13];
% grp2Retain = [1 2 3 6 7 12 13];
clc;
numNeighbors = 5;
grp2Retain = [8];
% Find clusters where the group has atleast 10% of cells
dMat = bsxfun(@rdivide,disMat,sum(disMat,1));
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
sorce2destination = nan(100,2);
% For each selected cluster find a cluster that is a direct link.
% If the direct link cluster is already selected, then ignore else
% add it to cls2Retain
cnt = 1;
for i = 1:numel(cls2Retain) 
    for j = i+1:numel(cls2Retain)
        sorce2destination(cnt,1) = cls2Retain(i);
        sorce2destination(cnt,2) = cls2Retain(j);
        cnt = cnt+1;
    end
end


% Find Connected clusters
for i = 1:numel(cls2Retain)
    ii = rowColumnClusterIndex(:,1) == cls2Retain(i);
    ii = or(ii,rowColumnClusterIndex(:,2) == cls2Retain(i));
    jj = rowColumnClusterIndex(:,3) >= 0;
    kk = and(ii,jj);
    indexValue = find(kk);
    cTmp = costValue(kk);
    [~,I] = sort(cTmp);
    indexValue = indexValue(I);
    if(numel(indexValue)>numNeighbors)
        kk(indexValue(numNeighbors+1:end)) = false;
    
    end
    
    
    if(sum(kk)>0)
        sorce2destination(cnt:cnt+sum(kk)-1,1) = repmat(cls2Retain(i),sum(kk),1);
        ix = rowColumnClusterIndex(kk,1:2);
        ix = ix(:);
        ix = ix(ix ~= cls2Retain(i));
        sorce2destination(cnt:cnt+sum(kk)-1,2) = ix;
        cnt = cnt+sum(kk);
    end    
end

ii = sum(isnan(sorce2destination),2)>0;
sorce2destination = sorce2destination(~ii,:);
ii = sorce2destination(:,1) ~= sorce2destination(:,2);
sorce2destination = sorce2destination(ii,:);
sorce2destination =  sort(sorce2destination,2);
% [~,I] = min(sorce2destination,[],2);
sorce2destination = unique(sorce2destination,'rows');
[~,I] = sort(sorce2destination(:,1));
sorce2destination = sorce2destination(I,:);
disp('###')

% Between Source and destination find the index of the clusters
clusterIndicestoRetain = nan(size(sorce2destination,1),1);
cnt = 1;
for i = 1:size(sorce2destination,1)
    ii = rowColumnClusterIndex(:,1) == sorce2destination(i,1);
    jj = rowColumnClusterIndex(:,2) == sorce2destination(i,2);
    kk = logical(ii.*jj);
    pIndex = cc{kk,1};
    m = intersect(pIndex(2:end-1),uIndx);
    if(~isempty(m))
        clusterIndicestoRetain(cnt:cnt+numel(m)-1,1) = m;
        cnt = cnt+numel(m);
    end
end
clusterIndicestoRetain = clusterIndicestoRetain(~isnan(clusterIndicestoRetain));
clusterIndicestoRetain= [clusterIndicestoRetain;sorce2destination(:)];
clusterIndicestoRetain = sort(unique(clusterIndicestoRetain));

% Get connection data based on clusters to retain
tmpPath = cell(10,1);
tmpIndex = nan(10,2);
for i = 1:size(sorce2destination,1)
    ii = rowColumnClusterIndex(:,1) == sorce2destination(i,1);
    jj = rowColumnClusterIndex(:,2) == sorce2destination(i,2);
    kk = logical(ii.*jj);
    pIndex = cc{kk,1};
    for k = 2:numel(pIndex)
        adjacencyMatrixGroup(pIndex(k-1),pIndex(k)) = 1;
        adjacencyMatrixGroup(pIndex(k),pIndex(k-1)) = 1;
    end
end
% tmpPath = tmpPath(1:cnt-1,1);
% tmpIndex = tmpIndex(1:cnt-1,:);
%

nClsIndex = false(numel(clusterIndex),1);
for i = 1:numel(clusterIndicestoRetain)
    ii = clusterIndex==clusterIndicestoRetain(i);
    nClsIndex(ii) = true;
end
nCindex = clusterIndex;
nCindex(~nClsIndex) = 0;
viewMSTPie2(redData,nCindex,mGrp,...
                                    maps,adjacencyMatrixGroup,edgeWidthMatrix,...
                                    controlNames);
%% Retain Groups and Plot MST between clusters without including intermediate
% Points
clc;
grp2Retain = [8];
% Find clusters where the group has atleast 10% of cells
dMat = bsxfun(@rdivide,disMat,sum(disMat,2));
dMat = dMat>.2;
dMat = dMat(:,grp2Retain);
cls2Retain = sum(dMat,2)>0;
adjacencyMatrixGroup = zeros(size(mCent,1));
edgeWeightMatrixGroup = inf(size(mCent,1));
edgeColor = zeros(size(mCent,1),size(mCent,1),3);
% grp2Retain = [4 5];
% Pull all clusters that have statistical enrichment
% pvalEnr = pval<=(.05/(numel(uIndx)*max(mGrp)));
% cls2Retain = sum(pvalEnr(grp2Retain,:))>0;
%  cls2Retain = sum(pvalEnr(grp2Retain,:))>0;
cls2Retain = sort(uIndx(cls2Retain));
sorce2destination = nan(100,2);
% For each selected cluster find a cluster that is a direct link.
% If the direct link cluster is already selected, then ignore else
% add it to cls2Retain
cnt = 1;
for i = 1:numel(cls2Retain) 
    for j = i+1:numel(cls2Retain)
        sorce2destination(cnt,1) = cls2Retain(i);
        sorce2destination(cnt,2) = cls2Retain(j);
        cnt = cnt+1;
    end
end


% Find Connected clusters
numNeighbors = 5;
for i = 1:numel(cls2Retain)
    ii = rowColumnClusterIndex(:,1) == cls2Retain(i);
    ii = or(ii,rowColumnClusterIndex(:,2) == cls2Retain(i));
    jj = rowColumnClusterIndex(:,3) == 0;
    kk = and(ii,jj);
    indexValue = find(kk);
    cTmp = costValue(kk);
    [~,I] = sort(cTmp);
   
    indexValue = indexValue(I);
    if(numel(indexValue)>numNeighbors)
        kk(indexValue(numNeighbors+1:end)) = false;
    
    end
    if(sum(kk)>0)
        sorce2destination(cnt:cnt+sum(kk)-1,1) = repmat(cls2Retain(i),sum(kk),1);
        ix = rowColumnClusterIndex(kk,1:2);
        ix = ix(:);
        ix = ix(ix ~= cls2Retain(i));
        sorce2destination(cnt:cnt+sum(kk)-1,2) = ix;
        cnt = cnt+sum(kk);
    end    
end

ii = sum(isnan(sorce2destination),2)>0;
sorce2destination = sorce2destination(~ii,:);
ii = sorce2destination(:,1) ~= sorce2destination(:,2);
sorce2destination = sorce2destination(ii,:);
sorce2destination =  sort(sorce2destination,2);
% [~,I] = min(sorce2destination,[],2);
sorce2destination = unique(sorce2destination,'rows');
[~,I] = sort(sorce2destination(:,1));
sorce2destination = sorce2destination(I,:);
disp('###')

% Between Source and destination find the index of the clusters
clusterIndicestoRetain = nan(size(sorce2destination,1),1);
cnt = 1;
for i = 1:size(sorce2destination,1)
    ii = rowColumnClusterIndex(:,1) == sorce2destination(i,1);
    jj = rowColumnClusterIndex(:,2) == sorce2destination(i,2);
    kk = logical(ii.*jj);
    pIndex = cc{kk,1};
    m = intersect(pIndex(2:end-1),uIndx);
    if(~isempty(m))
        clusterIndicestoRetain(cnt:cnt+numel(m)-1,1) = m;
        cnt = cnt+numel(m);
    end
end
clusterIndicestoRetain = clusterIndicestoRetain(~isnan(clusterIndicestoRetain));
clusterIndicestoRetain= [clusterIndicestoRetain;sorce2destination(:)];
clusterIndicestoRetain = sort(unique(clusterIndicestoRetain));
lineColor = zeros(size(sorce2destination,1),3);
% lineWidth = zeros(size(sorce2destination,1),3);
ccc = zeros(size(sorce2destination,1),1);
% Get connection data based on clusters to retain
tmpPath = cell(10,1);
tmpIndex = nan(10,2);
for i = 1:size(sorce2destination,1)
    ii = rowColumnClusterIndex(:,1) == sorce2destination(i,1);
    jj = rowColumnClusterIndex(:,2) == sorce2destination(i,2);
    kk = logical(ii.*jj);
    pIndex = cc{kk,1};
    ii = histc(mGrp(pIndex(2:end-1)),1:max(mGrp));
    mp = bsxfun(@times,maps,ii/sum(ii));
    mp = sum(mp);
    adjacencyMatrixGroup(sorce2destination(i,1),sorce2destination(i,2)) = 1;
    adjacencyMatrixGroup(sorce2destination(i,2),sorce2destination(i,1)) = 1;
    edgeWeightMatrixGroup(sorce2destination(i,1),sorce2destination(i,2)) = costValue(kk,1);
    edgeWeightMatrixGroup(sorce2destination(i,2),sorce2destination(i,1)) = costValue(kk,1);
    edgeColor(sorce2destination(i,1),sorce2destination(i,2),:) = reshape(mp,1,1,3);
    edgeColor(sorce2destination(i,2),sorce2destination(i,1),:) = reshape(mp,1,1,3);
    ccc(i) = costValue(kk,1);
end

nClsIndex = false(numel(clusterIndex),1);
for i = 1:numel(clusterIndicestoRetain)
    ii = clusterIndex==clusterIndicestoRetain(i);
    nClsIndex(ii) = true;
end
nCindex = clusterIndex;
nCindex(~nClsIndex) = 0;
viewMSTPieCluster2Cluster(redData,nCindex,mGrp,...
                                    maps,adjacencyMatrixGroup,edgeWeightMatrixGroup,...
                                    edgeColor,controlNames);


                                
                                
                                
                                
                                
                                
                                
                                
                                
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
