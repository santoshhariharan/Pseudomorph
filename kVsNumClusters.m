%% PSEUDOMORPH - CODE
% Read data
% Filter Data
%     - For Intensity
%     - For Morphology
%     - For Focus
%     - For expression
% Cluster data (Repeat)
%     - Per control K = 5
%     - Create MST for each run by jaccard
% Save multiple MST
%
clear;clc;close all;
warning('off','all');
%% Define Input Values
% 
% Module - 1:
% Define Variables 
featureReduction = false; % Feature reduction
numFeatures = 30;
pth='F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data';
load(fullfile(pth,'parameters.mat'));% Load parameter file
param.rootpath = pth;
intensityFeature = 'Ch2_INT_Cell_intensity';
intFeat = strcmpi(intensityFeature,param.datahdr);
% roiIntFeat = strcmpi('Ch2_MOR_cell_ROI_AvgIntensity',param.datahdr);
nucAreaFeat = strcmpi('Ch1_MOR_Nucleus_area',param.datahdr);
cellAreaFeat = strcmpi('Ch1_MOR_Cytoplasm_area',param.datahdr);
filePrefix = '.txt';
maxMinTType = false; % Normalization type
mxRw = 1000000;
numRandPrc = .7;
columnForControls = 9;
columnForOrganelle = 10;


fprintf('Starting pseudomorph\n');
fNames = dir(pth);
% featureReduction = true;
% clustersPerLandmark = true;


% Module 3: Read & Load data after filtering 
allD = zeros(mxRw,sum(param.datafeat));
allInten = zeros(mxRw,1);
allMorRatio = zeros(mxRw,1);
allTxt = cell(mxRw,1);
allTxtOrg = cell(mxRw,1);
% allMorIntensity = zeros(mxRw,1);
cnt = 1;
fprintf('Reading Data................');
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
    
    D.data = D.data(ii,:);
    D.textdata = D.textdata(ii,:);
    allInten(cnt:cnt+size(D.data,1)-1,:) = D.data(:,intFeat); 
%     allMorIntensity(cnt:cnt+size(D.data,1)-1,:) = D.data(:,roiIntFeat);
    allMorRatio(cnt:cnt+size(D.data,1)-1,:) = D.data(:,nucAreaFeat)./...
                                            (D.data(:,cellAreaFeat)+D.data(:,nucAreaFeat));    
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
%     allMorIntensity = allMorIntensity(1:cnt-1,:);
    allTxtOrg = allTxtOrg(1:cnt-1,:);
end

fprintf('\n');
% Remove nan entries
ii = sum(isnan(allD),2) ==0;
allD = allD(ii,:);
allInten = allInten(ii,:);
allTxt = allTxt(ii,:);
allMorRatio = allMorRatio(ii,:);
% allMorIntensity = allMorIntensity(ii,:);
allTxtOrg= allTxtOrg(ii,:);
fprintf('\nRemoved NAN %i\n',sum(~ii));

% Remove incorrectly segmented & low intensity Objects

ii = allMorRatio <= .5 & allMorRatio >= .2; % Value of this needs optimization
allD = allD(ii,:);
allTxt = allTxt(ii,:);
allTxtOrg = allTxtOrg(ii,:);
allInten = allInten(ii,:);
fprintf('#Cells Removed 4 Morphology %i of %i\n',sum(~ii),numel(ii));

% Normalization type
if(maxMinTType)
    minD = quantile(allD,.01);
    maxD = quantile(allD,.99);
    allD = bsxfun(@minus,allD, minD);
    allD = bsxfun(@rdivide,allD,maxD - minD);
    allInten = (allInten - min(allInten))./(max(allInten) - min(allInten));
else
    allInten = zscore(allInten);
    meanD = mean(allD);
    stdD = std(allD);
    allD = zscore(allD);    
end


% Remove intensity correlated features
rho = corr(allD,allInten);
newHeader = param.datahdr(1,param.datafeat);
ii = rho>-.5 & rho < .5;% Retain columns between -.5 & 0,5
allD = allD(:,ii);
fprintf('Features removed due to correlation\n');
fprintf('%s\n',newHeader{1,~ii});
newHeader = newHeader(1,ii);

% Remove Lower 1% and upper 1% data for each control
uControls = unique(allTxt(:,1));
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
% allMorIntensity = allMorIntensity(jj,:);
allTxtOrg = allTxtOrg(jj,:);
% param.meaninc = mean(allD);
% param.varinc = var(allD);
fprintf('#Cells removed by lower-upper quartile %i\n',sum(~jj));

% Remove features having 75% same data
[~,F] = mode(allD,1);
F= F./size(allD,1);
ii = F<.75;
newHeader = newHeader(1,ii);
allD = allD(:,ii);
% Feature Selection/Reduction
% newHeader = param.datahdr(1,param.datafeat);

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
allD = allD(:,redFeatures);
% Print number of cells per control
for i = 1:numel(uControls)
    ii = (strcmpi(uControls{i,:},allTxt(:,1)));
    fprintf('%s\t: %d\n',uControls{i,:},sum(ii));    
end


clear D focus cnt tok iFiles mxRw 
clear allMorRatio
clear ii jj kk rho mxRw 
clear intensityFeature intFeat roiIntFeat nucAreaFeat cellAreaFeat
clear fNames filePrefix randPrc allMorRatio allInten

%% Pick samples from each control randomly
% Test k values against number of clusters for each group using jaccard
% similarity
numK = [5:5:50];
numRandPrc = .7;
getEqualSamples  = true;
gps = getGroupIndices(allTxt,unique(allTxt));
numRpt = 1;
samIndex = false(numel(gps),1);
minSamplePerControl = 20000;
k = 5;
graphType = 'jaccard';
% Pick minimum set of samples
for i = 1:max(gps)
    minSamplePerControl = min(minSamplePerControl,sum(gps ==i ));
end
minSamplePerControl = floor(numRandPrc*minSamplePerControl);
fprintf('minimum samples per control %i\n',minSamplePerControl);
clsNum = zeros(max(gps),numel(numK));
for i = 1:max(gps)
    ii = find(gps == i);
    samIndex(ii(randperm(numel(ii),minSamplePerControl)),1) = true;
end

if(getEqualSamples)
    nGps = gps(samIndex);
    data4Clustering = allD(samIndex,:);
else
    nGps = gps;
end
maxCnt = max(nGps)*numel(numK);
h = waitbar(0,'Running K vs Cls calibration');
cnt = 0;
for iControl = 1:max(nGps)
    ii = nGps==iControl;
    fprintf('%s\n - %d\n',uControls{iControl,:},sum(ii));
    if(getEqualSamples)
        grpData = data4Clustering(ii,:);
    else
        grpData = allD(ii,:);
    end
    for j = 1:numel(numK)
        indx = phenograph(grpData,numK(j),'graphtype',graphType);
        clsNum(iControl,j) = numel(unique(indx));
        cnt = cnt+1;
        waitbar(cnt/maxCnt,h);
        
    end
end
disp('Done');
close(h);
clear grpData indx uIndx ii data4Clustering nGps;
clear iControl graphType samIndex
%%
