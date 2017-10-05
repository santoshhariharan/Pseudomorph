%%
% 
% Module - 1:
% Define Variables 
% Clear
clear;clc;close all;
% Define input variables
centroidMatFilename = 'centroidPerControl_Level2_RedFeature_5K_20K.mat';
load(centroidMatFilename);
oPFilename = 'clusterOP.txt';
writePerControl = true;
mCent = mCent2;
mGrp = mGrp2;
% Well-No. Plate-ID Image-No. Control Field-of-View X-Coord Y-Coord Classes
%%
fprintf('Starting pseudomorph\n');
pth='F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data';
load(fullfile(pth,'parameters.mat'));% Load parameter file
param.rootpath = pth;
intensityFeature = 'Ch2_INT_Cell_intensity';
intFeat = strcmpi(intensityFeature,param.datahdr);
nucAreaFeat = strcmpi('Ch1_MOR_Nucleus_area',param.datahdr);
cellAreaFeat = strcmpi('Ch1_MOR_Cytoplasm_area',param.datahdr);
filePrefix = '.txt';
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
% allTxtOrg = cell(mxRw,1);
txt4Writing = cell(mxRw,14);
cnt = 1;
fprintf('Completed Reading................');
for iFiles = 3:size(fNames,1)
    fprintf('\b\b\b\b\b\b\b\b\b%8.3f%%',iFiles*100./size(fNames,1));    
    if(fNames(iFiles).isdir)
        continue;
    end
    tok = regexpi(fNames(iFiles).name,'tacb5');
    if(~isempty(tok))
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
    
    ii = sum(isnan(D.data),2) == 0;
    
    D.data = D.data(ii,:);
    D.textdata = D.textdata(ii,:);
    
    allInten(cnt:cnt+size(D.data,1)-1,:) = D.data(:,intFeat); 
%     allMorIntensity(cnt:cnt+size(D.data,1)-1,:) = D.data(:,roiIntFeat);
    allMorRatio(cnt:cnt+size(D.data,1)-1,:) = D.data(:,nucAreaFeat)./...
                                            (D.data(:,cellAreaFeat)+D.data(:,nucAreaFeat));    
    allD(cnt:cnt+size(D.data,1)-1,:) = D.data(:,param.datafeat);
    allTxt(cnt:cnt+size(D.data,1)-1,:)= D.textdata(:,columnForControls);
    txt4Writing(cnt:cnt+size(D.data,1)-1,:) = D.textdata(:,1:14);
%     allTxt(cnt:cnt+size(D.data,1)-1,:)= D.textdata(:,:);
%     allTxtOrg(cnt:cnt+size(D.data,1)-1,:)= D.textdata(:,columnForOrganelle);
    cnt = cnt+size(D.data,1);
end
if(cnt<mxRw)
    allD = allD(1:cnt-1,:);
    allInten = allInten(1:cnt-1,:);
    allTxt = allTxt(1:cnt-1,:);
    allMorRatio = allMorRatio(1:cnt-1,:);
%     allMorIntensity = allMorIntensity(1:cnt-1,:);
%     allTxtOrg = allTxtOrg(1:cnt-1,:);
    txt4Writing = txt4Writing(1:cnt-1,:);
end

fprintf('\n');
% Remove nan entries
ii = sum(isnan(allD),2) ==0;
allD = allD(ii,:);
allInten = allInten(ii,:);
allTxt = allTxt(ii,:);
allMorRatio = allMorRatio(ii,:);
txt4Writing = txt4Writing(ii,:);
% allMorIntensity = allMorIntensity(ii,:);
% allTxtOrg= allTxtOrg(ii,:);
fprintf('\nRemoved NAN %i\n',sum(~ii));

ii = allInten >100;
allD = allD(ii,:);
allInten = allInten(ii,:);
allTxt = allTxt(ii,:);
allMorRatio = allMorRatio(ii,:);
txt4Writing = txt4Writing(ii,:);


% Remove incorrectly segmented & low intensity Objects
ii = allMorRatio <= .5 & allMorRatio >= .2; % Value of this needs optimization
allD = allD(ii,:);
allTxt = allTxt(ii,:);
% allTxtOrg = allTxtOrg(ii,:);
allInten = allInten(ii,:);
txt4Writing = txt4Writing(ii,:);
% allMorRatio = allMorRatio(ii,:);
% allMorIntensity = allMorIntensity(ii,:);
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
corrFeat = rho>-.5 & rho < .5;% Retain columns between -.5 & 0,5
allD = allD(:,corrFeat);
meanD = meanD(:,corrFeat);
stdD = stdD(:,corrFeat);
fprintf('#Features removed due to correlation %i\n',sum(~corrFeat));
fprintf('%s\n',newHeader{1,~corrFeat});
newHeader = newHeader(1,corrFeat);

% Remove Lower 1% and upper 1% data for each control
uControls = unique(allTxt);
jj = false(size(allTxt,1),1);
for i = 1:numel(uControls)
    ii = find(strcmpi(uControls{i,:},allTxt));
    kk = allInten(ii,1)>quantile(allInten(ii,1),.05) &...
                    allInten(ii,1)<quantile(allInten(ii,:),.95);
    jj(ii(kk)) = true;
end
allD = allD(jj,:);
allTxt = allTxt(jj,:);
% allTxtOrg = allTxtOrg(jj,:);
txt4Writing = txt4Writing(jj,:);
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
featureReduction = true;
numFeatures = 60;
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
%% Read centroid data
% load(centroidMatFilename)
numNeigh = 5;
allIdx = nan(numel(mGrp),numNeigh);
class = nan(numel(mGrp),numNeigh);
if(writePerControl)
    
    uC = unique(allTxt);
    grp = getGroupIndices(allTxt,uC);
    cnt = 1;
    for i = 1:max(grp)
        ii = find(grp == i);
        tmpD = allD(ii,:);
        idx = knnsearch(tmpD,mCent(mGrp==i,:),'K',numNeigh);
        allIdx(cnt:cnt+size(idx,1)-1,:) = ii(idx);
        class(cnt:cnt+size(idx,1)-1,:) = repmat([cnt:cnt+size(idx,1)-1]',1,numNeigh);
        cnt = cnt+size(idx,1);
%         lvl1 = knnclassify(tmpD, mCent(mGrp==i), 1:size(mCent,1));
    end
    
   ii =  sum(isnan(allIdx),2)==0;
   allIdx = allIdx(ii,:);
   class = class(ii,:);
    
else
    allIdx = knnsearch(allD,mCent,'K',numNeigh);
    class = repmat([1:numel(mGrp)]',1,numNeigh);
end
allIdx = allIdx(:);
class = class(:);
disp('####')
%%
writeClusters( class,txt4Writing(allIdx,:),oPFilename);



