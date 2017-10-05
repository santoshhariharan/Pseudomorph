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
%% Read data from files
% 
% Module - 1:
% Define Variables 

fprintf('Starting pseudomorph\n');
% Define Inputs:
pth='F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data';
featureReduction = true;
if(featureReduction)
    featStatus = 'True';
    numFeatures = 30;
else
    featStatus = 'False';
    numFeatures = 160;
end

columnForControls = 9;
columnForOrganelle = 10;
fprintf('Path: %s\n',pth);
fprintf('Feature reduction: %s\n',featStatus);
fprintf('---Number Feature: %i\n',numFeatures);
fprintf('Control column: %i\n',columnForControls);
load(fullfile(pth,'parameters.mat'));% Load parameter file
param.rootpath = pth;
intensityFeature = 'Ch2_INT_Cell_intensity';
intFeat = strcmpi(intensityFeature,param.datahdr);
nucAreaFeat = strcmpi('Ch1_MOR_Nucleus_area',param.datahdr);
cellAreaFeat = strcmpi('Ch1_MOR_Cytoplasm_area',param.datahdr);
filePrefix = '.txt';
fNames = dir(pth);


% Module 3: Read & Load data after filtering
fprintf('Module 3.......\n');
mxRw = 1000000; 
allD = zeros(mxRw,sum(param.datafeat));
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
%     allTxt(cnt:cnt+size(D.data,1)-1,:)= D.textdata(:,:);
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
% Remove Artifacts and noise

% Remove NaN entries
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
% allMorRatio = allMorRatio(ii,:);
% allMorIntensity = allMorIntensity(ii,:);
fprintf('#Cells Removed 4 Morphology %i of %i\n',sum(~ii),numel(ii));


% Remove cells with low intensity
ii = allInten > 100;
allD = allD(ii,:);
allTxt = allTxt(ii,:);
allTxtOrg = allTxtOrg(ii,:);
allInten = allInten(ii,:);
fprintf('#Cells Removed 4 Intensity %i of %\n',sum(~ii));

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
allTxtOrg = allTxtOrg(jj,:);
allInten = allInten(jj,:);
% param.meaninc = mean(allD);
% param.varinc = var(allD);
fprintf('#Cells removed by lower-upper quartile %i\n',sum(~jj));


% Remove intensity correlated features
rho = corr(allD,allInten);
newHeader = param.datahdr(1,param.datafeat);
ii = find(param.datafeat);
corrFeat = rho>-.5 & rho < .5;% Retain columns between -.5 & 0,5
param.datafeat(1,ii(~corrFeat)) = false;
allD = allD(:,corrFeat);
fprintf('#Features removed due to correlation %i\n',sum(~corrFeat));
fprintf('%s\n',newHeader{1,~corrFeat});
newHeader = newHeader(1,corrFeat);


% Remove features having 75% same data
[~,F] = mode(allD,1);
F= F./size(allD,1);
jj = F<.75;
ii = find(param.datafeat);
param.datafeat(1,ii(~jj)) = false;

newHeader = newHeader(1,jj);
allD = allD(:,jj);
% Feature Selection/Reduction
% newHeader = param.datahdr(1,param.datafeat);


if(numel(newHeader)<=numFeatures)
    featureReduction = false;
    fprintf('TURNED OFF FEATURE REDUCTION\n');
end
if(featureReduction)    
    redFeatures = unsupervisedGreedyFS(allD,numFeatures);
else
    redFeatures = true(1,size(allD,2));
end

ii = find(param.datafeat);
param.datafeat(1,ii(~redFeatures)) = false;
fprintf('Features Chosen\n')
fprintf('%s\n',newHeader{1,redFeatures});
% [cf]= princomp(allD(:,redFeatures));
allD = allD(:,redFeatures);
% Print number of cells per control
for i = 1:numel(uControls)
    ii = (strcmpi(uControls{i,:},allTxt));
    fprintf('%s\t: %d\n',uControls{i,:},sum(ii));    
end

% Normalization
meanD = mean(allD);
stdD = std(allD);
% allD = zscore(allD);  
minD = min(allD);
maxD = max(allD);
allD = bsxfun(@minus,allD,minD);
allD = bsxfun(@rdivide,allD,maxD-minD);

clear D focus cnt tok iFiles mxRw 
clear allMorRatio allInten
clear ii jj kk rho mxRw 
clear intensityFeature roiIntFeat
clear fNames filePrefix randPrc



%% Select Samples for Vieweing data
vData = [];
vGrp = [];
mc = [];
grp = getGroupIndices(allTxt,unique(allTxt));
for i =1:max(grp)
    ii = find(grp == i);
    idx = knnsearch(allD(ii,:),mean(allD(ii,:)),'k',100);
    idx = idx';
    vData = [vData;allD(ii(idx),:)];
    vGrp = [vGrp;grp(ii(idx,1))];
    mc = [mc;mean(allD(ii,:))];
end
pp= pdist2(mc,mc);
%% tSNE
redData = compute_mapping(vData,'Sammon',3);
disp('$$$$')
%%
figure;hold on;
for i = 1:max(grp)
    ii = vGrp == i;
%     plot(redData(ii,1),redData(ii,2),'o','MarkerFacecolor',map(i,:),...
%         'MarkerEdgeColor','None');
    plot3(redData(ii,1),redData(ii,2),redData(ii,3),'o','MarkerFacecolor',map(i,:),...
        'MarkerEdgeColor','None');
end
hold off;
legend(unique(allTxt));

