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

% Module - 1:
% Define Variables 
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
%% Feature Selection/Reduction
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
options.kcls = 10;
options.k = 5;
options.sampling = 'random';
options.repeat = 1;
gps = getGroupIndices(allTxt,unique(allTxt));
uG = unique(gps);
mCent = nan(500,size(allD,2));
mGrp = nan(500,1);
numPointsPerCluster = nan(500,1);
grpCentroids = zeros(numel(uG),size(allD,2));

cnt = 1;
for i = 1:numel(uG)
    ii = gps == uG(i);
    grpData = allD(ii,:);
    options.grouped = ones(sum(ii),1);
    grpCentroids(i,:) = mean(grpData);
    [ indx, mIndx, objectAssign ] = consensusClustering( grpData,options);
    uInd = unique(indx);
    for j = 1:numel(uInd)
       mCent(cnt,:)= mean(grpData(indx == uInd(j),:));
       mGrp(cnt,1) = i;
       numPointsPerCluster(cnt,1) = sum(indx == uInd(j));
       cnt = cnt+1;
    end
end
ii = ~isnan(mGrp);
mCent = mCent(ii,:);
mGrp = mGrp(ii,1);
numPointsPerCluster = numPointsPerCluster(ii,1);
clear options grpData i j uInd indx cnt ii;
%% Number of Points per cluster
uIndx = unique(indx);
for i = 1:numel(uIndx)
    fprintf('Cluster%i\t%.0f\n',i,sum(indx == uIndx(i)));
end


%% Compute clustering of centroids
indxN = phenograph(mCent,5);
%% Compute t-SNE

redPlot = compute_mapping(mCent,'t-SNE',2);

%% Plot data
uGrp = unique(mGrp);
figure;hold on;
for i = 1:numel(uGrp)
    plot(redPlot(mGrp == uGrp(i),1),redPlot(mGrp == uGrp(i),2),'o',...
                'Markerfacecolor',param.maps(i,:),'Markeredgecolor','none');
end
hold  off;
legend(uControls);
%%

uIndx = unique(indxN);
mps= jet(numel(uIndx));
figure;hold on;
for i = 1:numel(uIndx)
    plot(redPlot(indxN == uIndx(i),1),redPlot(indxN == uIndx(i),2),'o',...
                'Markerfacecolor',mps(i,:),'Markeredgecolor','none');
end
hold  off;
% legend(uControls);

%% Compute weighted distance between points

% Compute euclidean distance between groups
weights = pdist2(grpCentroids,grpCentroids);
wTmp = weights(:);
wTmp = wTmp(wTmp>0);
sigma = std(wTmp);
weights = weights.^2;
weights = weights/sigma.^2;
weights = exp(weights);


[m,n] = size(mCent);
centSim = nan(m*(m-1)/2,3);cnt = 1;
for i = 1:m
    for j = i+1:m
        centSim(cnt,1) = j;centSim(cnt,2) = i;
        centsim(cnt,3) = pdist2(mCent(i,:),mCent(j,:));       
        centsim(cnt,3) = centsim(cnt,3).*weights(mGrp(1,1),mGrp(j,1));
        cnt = cnt+1;        
    end
end
        
[~,~,T] = kruskal(centSim(:,[1 2]),centSim(:,3));

% Compute Path length

% Cluster based on path length
    



%% Compute Centroids from level 2
indx2 = phenograph( clsCentroidLevel1, 3,'distance','Jaccard');
uIndx = unique(indx2);
clsCentroidLevel2 = zeros(numel(uIndx),size(clsCentroidLevel1,2));
for i = 1:numel(uIndx)
    ii = indx2 == uIndx(i);
    clsCentroidLevel2(i,:) = mean(clsCentroidLevel1(ii,:));
end


%% Compute MST @ different levels
% Main Level

D = pdist2(controlCentroids,controlCentroids);
[~,~,xst] = kruskal(1-eye(size(D,1)),D);
writeCytoscapeNetworkFile([xst D(xst(:,1),xst(:,2))],...
    fullfile('F:\Projects\Proteinlocalization\PseudoMorph','Level1Network') );







%%

% [cf] = princomp(allD);
% cf = cf(:,1:2);
% scr = clsCentroidLevel1*cf;

scr=compute_mapping(clsCentroidLevel1, 't-SNE', 2);
%% Plot Filter with less than 50

orderControls = {'CB5-ER';'ER-PRO';'TABIK';'ERGIC';'GOLGI-GT';'GOLGIN';'MAO';'MITO-CCO';...
    'METAXIN';'TACB5';'PQC-PSS1';'RAB5A';'RAB7A';'FL_VAMP2';'FL-VAMP5';...
    'CYTO';'LAMINA';'PEROXISOME-1';'LAMP1'};
orderControls = unique(orderControls);
ii = elementsPerCluster >=50;
mSize = 100*ones(sum(ii),1);        
clsCentroidLevel1 = clsCentroidLevel1(ii,:);
controlCategory = controlCategory(ii,:);
elementsPerCluster = elementsPerCluster(ii,:);
mSize = ((.8)*((elementsPerCluster - min(elementsPerCluster))./...
            (max(elementsPerCluster) - min(elementsPerCluster))) + .1).*mSize;
        
% scr=compute_mapping(clsCentroidLevel1, 'PCA', 2);
kk = false(numel(orderControls),1);
figure;hold on;
for i = 1:numel(uControls)
    ii = strcmpi(controlCategory,uControls{i,:});
    jj = strcmpi(orderControls,uControls{i,:});
    kk(jj) = true;
    %     plot3(scr(ii,1),scr(ii,2),scr(ii,3),'o','MarkerFaceColor',param.maps(i,:),...
    %         'MarkerEdgeColor','none');
    scatter(scr(ii,1),scr(ii,2),mSize(ii,1),param.maps(jj,:),...
                    'filled');
end
hold off;
legend(uControls);title('Level 1 Centroids');
clear i iCnt ii iIndx tmp
%% Assign all data to centroids

allIndex = knnsearch(clsCentroidLevel1,allD);
%% Load Bio Mapping
load('biologicalMapping.mat');
biologicalMap = logical(biologicalMap);
bMap = inf*ones(size(clsCentroidLevel1,1),size(clsCentroidLevel1,1));
for i = 1:size(biologicalMap,1)
    ii = strcmpi(controlCategory,bioMapControlNames{i,1});
    
    mp2 = bioMapControlNames(biologicalMap(i,:));
    for j = 1:numel(mp2)
        if(isempty(mp2{j}))
            continue;
        end
        jj = (strcmpi(controlCategory,mp2{j}));
        bMap(ii,jj) = 0;
    end
end
disp('@@@@')
%% Find minimum spanning tree for each landmark
clc;
uC = unique(allTxt);
matA = nan(size(clsCentroidLevel1,1),1);% Number of edges stored
neighbourMatrix = nan(2000,2);
cnt = 1;nn = 1;
figure; hold on;
for i = 1:numel(uC)
    ii = find(strcmpi(uC{i,:},controlCategory));
    tmp = scr(ii,:);
    jj = strcmpi(orderControls,uC{i,:});
    D = pdist2(scr(ii,:),scr(ii,:),'euclidean');
    scatter(scr(ii,1),scr(ii,2),mSize(ii,1),param.maps(jj,:),...
                    'filled');
    [~,st,adj] = kruskal(1-eye(numel(ii)),D);
    for k = 1:size(st,1)
        coord = tmp(st(k,:),:)';
        line(coord(1,:),coord(2,:),'Linewidth',.2,'Color',[.4 .4 .4]);
    end
    matA(cnt:cnt+numel(ii)-1,1) = sum(adj,2);
    neighbourMatrix(nn:nn+size(st,1)-1,1) = ii(st(:,1));
    neighbourMatrix(nn:nn+size(st,1)-1,2) = ii(st(:,2));
    cnt = cnt+numel(ii);nn = nn + size(st,1);
    
end
hold off;
% Create Map suxh that only free points are avilable to form edges
aa = inf*ones(size(clsCentroidLevel1,1),size(clsCentroidLevel1,1));
ii = matA == 1;
aa(ii,ii) = 0;
aa = diag(inf.*ones(size(aa,1),1),0)+aa;



% scr = compute_mapping(clsCentroidLevel1,'t-SNE',2);
%% Compute distribution
uC = unique(allTxt);
uIndxL2 = unique(allIndex);
cDistribution = zeros(numel(uC),numel(uIndxL2));
for j = 1:numel(uC)
    jj = strcmpi(allTxt,uC{j,:});
    for i = 1:numel(uIndxL2)
        ii = allIndex==uIndxL2(i);        
        cDistribution(j,i) = sum(jj.*ii);
    end
end
%
cDistribution = cDistribution';
cDistributionP = bsxfun(@rdivide,cDistribution,sum(cDistribution,1));
D = pdist2(cDistributionP,cDistributionP,'correlation');

%
D = (D - min(D(:)))./(max(D(:))-min(D(:)));
D1 = pdist2(scr,scr,'euclidean');
D1 = (D1 - min(D1(:)))./(max(D1(:))-min(D1(:)));
if(~exist('bMap','var'))
    [w,xst] = kruskal(1-eye(size(scr,1)),D1);
else
    mt = bMap+(D1)+D;
    [w,xst] = kruskal(1-eye(size(scr,1)),mt);
end
neighbourMatrix(nn:nn+size(xst,1)-1,:) = xst;
neighbourMatrix = neighbourMatrix(sum(~isnan(neighbourMatrix),2)==2,:);
h=viewMSTPie2(scr,...
    cDistributionP,param.maps(kk,:),uC,xst);
return;
