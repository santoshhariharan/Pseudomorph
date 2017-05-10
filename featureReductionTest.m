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
% allInten = allInten(jj,:);
% allMorRatio = allMorRatio(jj,:);
allMorIntensity = allMorIntensity(jj,:);
allTxtOrg = allTxtOrg(jj,:);
param.meaninc = mean(allD);
param.varinc = var(allD);

% Inputs
k = 10;
distanceTypeForSampling = 'Euclidean';
graphType = 'jaccard';
gps = getGroupIndices(allTxt,unique(allTxt));
uG = unique(gps);
% Print number of cells per control
for i = 1:numel(uControls)
    ii = (strcmpi(uControls{i,:},allTxt));
    fprintf('%s\t: %d\n',uControls{i,:},sum(ii));    
end
% Feature Selection/Reduction
% if(featureReduction)    


clear D focus cnt tok iFiles mxRw 
clear allMorRatio allInten
clear ii jj kk rho mxRw 
clear intensityFeature intFeat roiIntFeat nucAreaFeat cellAreaFeat
clear fNames filePrefix randPrc
% Perform clustering high number of centroids -
% Uneven number of samples
% Feature Selection/Reduction

numFeatures = [10:10:100];

allCls = zeros(size(allD,1),numel(numFeatures));
for jFeatures = 1:numel(numFeatures)
    fprintf('#Features %d\n',numFeatures(jFeatures));
    redFeatures = unsupervisedGreedyFS(allD,numFeatures(jFeatures));
    nAllD = allD(:,redFeatures);
    cnt = 1;
    mCent = nan(500,size(nAllD,2));
    
    
    for i = 1:numel(uG)
        ii = gps ==uG(i);
        grpData = nAllD(ii,:);
        [ indx ] = phenograph(grpData,k,'graphtype',graphType);
        uInd = unique(indx);
        for j = 1:numel(uInd)
            kk =  indx == uInd(j);
            mCent(cnt,:)= mean(grpData(kk,:));
            cnt = cnt+1;
        end
    end
    ii = sum(isnan(mCent),2)==0;
    mCent = mCent(ii,:);
    allCls(:,jFeatures) = knnclassify(nAllD,mCent,1:size(mCent,1));
end
clear options grpData i j uInd indx cnt ii;
clear comm kk grpInten
%% Compute Adjusted Rand index
clc;
aRandIndex = nan(100,3);
m = numel(numFeatures);
cnt = 1;
totalN = m.*(m-1)./2;
fprintf('Computing ARI.............');
for i = 1:m
    fprintf('\b\b\b\b\b\b\b\b%7.3f%%',cnt*100/totalN);
    for j = i+1:m
        if(i==j)
            continue;
        end
        aRandIndex(cnt,1) = numFeatures(i);
        aRandIndex(cnt,2) = numFeatures(j);
        aRandIndex(cnt,3) = adjRandIndex(allCls(:,i),allCls(:,j));
        cnt = cnt+1;
    end
    
end
fprintf('Done\n')
aRandIndex = aRandIndex(sum(isnan(aRandIndex),2)==0,:);

