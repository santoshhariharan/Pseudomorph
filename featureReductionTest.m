clear;clc;close all;
warning('off','all');
% numDimsVis = 2;
fprintf('Starting pseudomorph\n');
pth='F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data\Mito-ERData';
load(fullfile(pth,'parameters.mat'));% Load parameter file
param.rootpath = pth;
intensityFeature = 'Ch2_INT_Cell_intensity';
intFeat = strcmpi(intensityFeature,param.datahdr);
roiIntFeat = strcmpi('Ch2_MOR_cell_ROI_AvgIntensity',param.datahdr);
nucAreaFeat = strcmpi('Ch1_MOR_Nucleus_area',param.datahdr);
cellAreaFeat = strcmpi('Ch1_MOR_Cytoplasm_area',param.datahdr);
filePrefix = '.txt';
% randPrc = .1;
% numFeatures = 30;
fNames = dir(pth);
columnForControls = 9;
columnForOrganelle = 10;
% featureReduction = true;
% clustersPerLandmark = true;

maxMinTType = true;
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
    
    D.data = D.data(ii,:);
    D.textdata = D.textdata(ii,:);
    allInten(cnt:cnt+size(D.data,1)-1,:) = D.data(:,intFeat); 
    allMorIntensity(cnt:cnt+size(D.data,1)-1,:) = D.data(:,roiIntFeat);
%     pp = D.data(:,nucAreaFeat)./(D.data(:,cellAreaFeat)+D.data(:,nucAreaFeat));
%     if(sum(pp>10))
%         fprintf('%s -- %f\n',fNames(iFiles).name,sum(pp>10)/numel(pp));
%     end
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
fprintf('\nRemoved NAN\n');
% hist(allMorRatio,sqrt(size(allMorRatio,1)));
% re
% Remove incorrectly segmented & low intensity Objects
ii = allMorRatio <= .5;
allD = allD(ii,:);
allTxt = allTxt(ii,:);
allTxtOrg = allTxtOrg(ii,:);
allInten = allInten(ii,:);
% allMorRatio = allMorRatio(ii,:);
% allMorIntensity = allMorIntensity(ii,:);
fprintf('#Cells Removed 4 intensity %i of %i\n',sum(~ii),numel(ii));

% Normalization type
if(maxMinTType)
    minD = quantile(allD,.01,1);
    maxD = quantile(allD,.99,1);
    allD = bsxfun(@minus,allD, minD);
    allD = bsxfun(@rdivide,allD,maxD - minD);
    allInten = (allInten - min(allInten))./(max(allInten) - min(allInten));
else
    allInten = zscore(allInten);
    allD = zscore(allD);
end


% Remove intensity correlated features
rho = corr(allD,allInten);
newHeader = param.datahdr(1,param.datafeat);
ii = rho>-.5 & rho < .5;% Retain columns between -.5 & 0,5
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
%% Pick samples from each control randomly
gps = getGroupIndices(allTxt,unique(allTxt));
samIndex = false(numel(gps),1);
minSamplePerControl = 20000;
% Pick minimum set of samples
for i = 1:max(gps)
    minSamplePerControl = min(minSamplePerControl,sum(gps ==i ));
end
minSamplePerControl = floor(.3*minSamplePerControl);
for i = 1:max(gps)
    ii = find(gps == i);    
    samIndex(ii(randperm(numel(ii),minSamplePerControl))) = true;
end
fprintf('minimum samples per control %i\n',minSamplePerControl);

sampleData = allD(samIndex,:);
%% Perform clustering high number of centroids -
% Uneven number of samples
% Feature Selection/Reduction
k = 40;
graphType = 'Jaccard';
numFeatures = [10:10:160];
allCls = zeros(sum(samIndex),numel(numFeatures));
for jFeatures = 1:numel(numFeatures)
    fprintf('#Features %d\n',numFeatures(jFeatures));
    redFeatures = unsupervisedGreedyFS(sampleData,numFeatures(jFeatures));
    nAllD = sampleData(:,redFeatures);
    allCls(:,jFeatures) = phenograph(nAllD,k,'graphtype',graphType);
%     C = clsIn(nAllD);
%     pref = C.pmed-((C.pmed-C.pmin)/(2^6));   
%     allCls(:,jFeatures) = apcluster(C.S,pref);
end
clc;disp('DONE');
clear options grpData i j uInd indx cnt ii;
clear comm kk grpInten jFeatures redFeatures
%% Compute Adjusted Rand index
% clc;
% aRandIndex = nan(100,3);
aRandIndex = zeros(numel(numFeatures),numel(numFeatures));
m = numel(numFeatures);
% cnt = 1;
% totalN = m.*(m-1)./2;
fprintf('Computing ARI.............\n');
for i = 1:m
    for j = i+1:m
        aRandIndex(i,j) = adjRandIndex(allCls(:,i),allCls(:,j));
    end
end
aRandIndex = aRandIndex+aRandIndex';
aRandIndex = aRandIndex/2;
xx = diag(aRandIndex,-1);

