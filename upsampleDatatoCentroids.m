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
%%
% 
% Module - 1:
% Define Variables 

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
featureReduction = false;
numFeatures = 30;
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
clear allMorRatio
clear ii jj kk rho mxRw 
clear intensityFeature intFeat roiIntFeat nucAreaFeat cellAreaFeat
clear fNames filePrefix randPrc

%% Load centroid data
load('Centroid5K.mat');
% Upsample data
lvl1 = knnclassify(allD, mCent, 1:size(mCent,1));
%% Discard Centroids with less than 5% of data points
uL1 = unique(lvl1);
m = numel(lvl1);
fprintf('Cluster# \t Group \t #Percent\n');
fprintf('--------------------\n');
cnt = 0;
clusterPercent = zeros(numel(uL1),1);
for i = 1:numel(uL1)
    fprintf('%2.0f \t %2.0f \t %6.3f\n ',i,mGrp(uL1(i)),...
                        sum(lvl1==uL1(i))*100/m);
    clusterPercent(i) = sum(lvl1==uL1(i))*100/m;                
    cnt = cnt+sum(lvl1==uL1(i));
end
figure;hist(clusterPercent,sqrt(numel(clusterPercent)));
xlabel('% of Total Data');
%%
numDimsVis = 2;
pcaRedData = compute_mapping(mCent,'SNE',numDimsVis);
disp('@@@@@@@@@@@@');
%% Do leave one out & compute MST

allMST = zeros(size(mCent,1));
for i = 1:max(mGrp)
    mp = true(max(mGrp),1);
    mp(i) = false;
    ii = find(mGrp ~= i);
    tmpCent = pcaRedData(ii,:);
%     [idx, dis] = knnsearch(tmpCent,tmpCent,'K',kVal+1);
%     idx = idx(:,2:end);
%     dis = dis(:,2:end);
%     adjMat = zeros(size(tmpCent,1));
%     distanceMatrix = inf(size(tmpCent,1));
%     for j = 1:size(tmpCent,1)
%         adjMat(j,idx(j,:)) = 1;
%         distanceMatrix(j,idx(j,:)) = dis(j,:);
%     end
    adjMat = 1-eye(size(tmpCent,1));
    distanceMatrix = pdist2(tmpCent,tmpCent);
    [~,xst] = kruskal(adjMat, distanceMatrix);
    plotMSTFigure( pcaRedData(ii,:),xst,mGrp(ii,1),map(mp,:) );
    for j = 1:size(xst,1)
        allMST(ii(xst(j,1)),ii(xst(j,2))) = allMST(ii(xst(j,1)),ii(xst(j,2)))+1;
    end
end
[r,c] = find(allMST>0);
plotMSTFigure( pcaRedData,[r c],mGrp,map );
disp('Done');
% [idx, dis] = knnsearch(mCent,mCent,'K',kVal+1);





%% View data
numDimsVis = 2;
pcaRedData = compute_mapping(mCent,'t-SNE',numDimsVis);
disp('@@@@@@@@@@@@');
%
kVal = 5;
[idx, dis] = knnsearch(mCent,mCent,'K',kVal+1);
idx = idx(:,2:end);
[m, n] = size(mCent);
% Plot data
figure; hold on;
for iControl = 1:m
    line([pcaRedData(iControl,1) pcaRedData(idx(iControl,:),1)'],...
            [pcaRedData(iControl,2) pcaRedData(idx(iControl,:),2)'],...
            'Color',[.8 .8 .8]);
end

for iControl = 1:max(mGrp)
%     if(sum(rmGrp == iControl)>0)
%         continue;
%     end
    
    ii = mGrp == iControl;
    plot(pcaRedData(ii,1),pcaRedData(ii,2),'o','MarkerFaceColor',map(iControl,:),...
        'MarkerEdgeColor','None','MarkerSize',6);      
end
% legend(uControls);
% set(gca,'XTick',[]);set(gca,'YTick',[]);
hold off;
%% Visualize using Pie chart

%% 3D Plots
% numDimsVis = 2;
pcaRedData = compute_mapping(mCent,'t-SNE',3);
disp('@@@@@@@@@@@@');
%
kVal = 1;
[idx, dis] = knnsearch(mCent,mCent,'K',kVal+1);
idx = idx(:,2:end);
[m, n] = size(mCent);
% Plot data
figure; hold on;
for iControl = 1:m
    line([pcaRedData(iControl,1) pcaRedData(idx(iControl,:),1)'],...
            [pcaRedData(iControl,2) pcaRedData(idx(iControl,:),2)'],...
            [pcaRedData(iControl,3) pcaRedData(idx(iControl,:),3)'],...
            'Color',[.8 .8 .8]);
end

for iControl = 1:max(mGrp)
%     if(sum(rmGrp == iControl)>0)
%         continue;
%     end
    
    ii = mGrp == iControl;
    plot3(pcaRedData(ii,1),pcaRedData(ii,2),pcaRedData(ii,3),'o','MarkerFaceColor',map(iControl,:),...
        'MarkerEdgeColor','None','MarkerSize',6);      
end
% legend(uControls);
% set(gca,'XTick',[]);set(gca,'YTick',[]);
hold off;

