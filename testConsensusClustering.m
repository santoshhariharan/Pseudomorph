% Test Consensus Clustering 
clear;clc;
fprintf('Starting pseudomorph\n');
pth='F:\Projects\Proteinlocalization\PseudoMorph\TestDataForPseudoMorph';
load(fullfile(pth,'parameters.mat'));% Load parameter file
param.rootpath = pth;
intensityFeature = 'Ch2_INT_Cell_intensity';
intFeat = strcmpi(intensityFeature,param.datahdr);
filePrefix = 'Controls';
randPrc = .1;
% numFeatures = 30;
fNames = dir(pth);
columnForControls = 9;
% featureReduction = true;
clustersPerLandmark = true;


%% Module 3: Read & Load data after filtering
fprintf('Module 3.......\n');
mxRw = 1000000;
allD = zeros(mxRw,sum(param.datafeat));
allInten = zeros(mxRw,1);
allTxt = cell(mxRw,1);
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
    allD(cnt:cnt+size(D.data,1)-1,:) = D.data(:,param.datafeat);
    allTxt(cnt:cnt+size(D.data,1)-1,:)= D.textdata(:,columnForControls);
    cnt = cnt+size(D.data,1);
end
if(cnt<mxRw)
    allD = allD(1:cnt-1,:);
    allInten = allInten(1:cnt-1,:);
    allTxt = allTxt(1:cnt-1,:);
end
fprintf('\n');
ii = allInten>100;
allD = allD(ii,:);
allTxt = allTxt(ii,:);
allInten = allInten(ii,:);
% allInten = allInten(ii,:);
% allD = bsxfun(@minus,allD, param.meaninc(1,param.datafeat));
% allD = bsxfun(@rdivide,allD,sqrt(param.varinc(1,param.datafeat)));

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
param.meaninc = mean(allD);
param.varinc = var(allD);
allD = zscore(allD);
% allInten = zscore(allInten);
clear D focus cnt tok iFiles mxRw allInten
% Feature Selection/Reduction
newHeader = param.datahdr(1,param.datafeat);
featureReduction = false;
numFeatures = 20;
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


%% Perform consensus clustering
options.distance = 'Jaccard';
options.algorithm = 'phenograph';
options.pref = 'pmin6';
options.kcls = 50;
options.k = 5;
options.sampling = 'random';
options.repeat = 10;
groupNames = unique(allTxt);
gps = getGroupIndices(allTxt,groupNames);
groupNumber = 1;
kVal = [5:5:50];
clustersByConsensus = zeros(numel(kVal),1);
clustersByPhenographMean = zeros(numel(kVal),1);
clustersByPhenographDev = zeros(numel(kVal),1);
fprintf('Group : %s\n',groupNames{groupNumber,:});
fprintf('_______________________________________\n');
fprintf('Completed........................');
for i = 1:numel(kVal)
    fprintf('\b\b\b\b\b\b\b\b%7.3f%%',i*100/numel(kVal));
    ii = gps == groupNumber;
    grpData = allD(ii,:);
    options.grouped = ones(sum(ii),1);
    options.k = kVal(i);
    [ indx, mIndx, ~ ] = consensusClustering( grpData,options);
    clustersByConsensus(i) = numel(unique(indx));
    tmp = zeros(options.repeat,1);
    for j = 1:options.repeat
        tmp(j) = numel(unique(mIndx(:,j)));
    end
    clustersByPhenographMean(i) = mean(tmp);
    clustersByPhenographDev(i) = std(tmp);    
end
fprintf('\n');
clear options grpData i j uInd indx cnt ii;
clear kVal tmp indx mIndx objectAssign;