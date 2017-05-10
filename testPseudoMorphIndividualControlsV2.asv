% Pseudomorph: Per control
% Run pseudo morph on individual controls. Over cluster the data so we have
% more clusters than necessary. This is controlled by variable max cluster
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
fprintf('Starting pseudomorph')
clear;clc;
pth='F:\Projects\Proteinlocalization\PseudoMorph\TestDataForPseudoMorph';
load(fullfile(pth,'parameters.mat'));% Load parameter file
param.rootpath = pth;
cntrl = {'CB5-ER';'ER-PRO';'ERGIC';'FL-VAMP5';'FL-VAMP2';...
    'GOLGI-GT';'GOLGIN';'METAXIN';...
    'PQC-PSS1';'RAB5A';'RAB7A';'MAO';'MITO-CCO';...
    'TABIK'};
intensityFeature = 'Ch2_INT_Cell_intensity';
intFeat = strcmpi(intensityFeature,param.datahdr);
filePrefix = 'Controls';
randPrc = .1;
numFeatures = 30;
fNames = dir(pth);
columnForControls = 9;
%% Module 2 - Compute mean and variance
fprintf('Module 2.......\n');
md = computeMeanVarianceCls(pth,param);
param.meaninc = md.meaninc;
param.varinc = md.varinc;
param.maxColVal = md.maxColVal;
param.minColVal = md.minColVal;
clear md;clc;
fprintf('Completed Mean & Variance Computation\n');


%% Module 3: Read some Sample data @ random
fprintf('Module 3.......\n');
mxRw = 20000;
allD = zeros(mxRw,sum(param.datafeat));
allInten = zeros(mxRw,1);
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
    p = randperm(size(D.data,1),floor(randPrc.*size(D.data,1)));
    D.data = D.data(p,:);
    allInten(cnt:cnt+size(D.data,1)-1,:) = D.data(:,intFeat);
    allD(cnt:cnt+size(D.data,1)-1,:) = D.data(:,param.datafeat);
    cnt = cnt+size(D.data,1);
end
if(cnt<mxRw)
    allD = allD(1:cnt-1,:);
    allInten = allInten(1:cnt-1,:);
end
fprintf('\n')
allD = allD(allInten>100,:);
allD = bsxfun(@minus,allD, param.meaninc(1,param.datafeat));
allD = bsxfun(@rdivide,allD,sqrt(param.varinc(1,param.datafeat)));
redFeatures = unsupervisedGreedyFS(allD,numFeatures);
% [cf]= princomp(allD(:,redFeatures));
clear D focus cnt tok iFiles allD mxRw allInten



%% Module 4:
fprintf('Module 4.......\n');
maxClsPerControl = 100;distanceType = 'Euclidean';
fileMap=mapFilesToControls(pth,cntrl);
uControls = unique(fileMap(:,2));
bCnt = 1;clsCentroidLevel1 = zeros(maxClsPerControl.*numel(cntrl),sum(redFeatures));
controlCategory = cell(maxClsPerControl.*numel(cntrl),1);
for iCnt = 1:numel(uControls)
    fprintf('%s\n',uControls{iCnt,:});
    filesPerControl = fileMap(strcmpi(fileMap(:,2),uControls{iCnt,:}),1);
    mxRw = 50000;cnt = 1;
    allD = zeros(mxRw,sum(param.datafeat));
    allInten = zeros(mxRw,1);
    filesPerControl = filesPerControl(randperm(numel(filesPerControl)),:);    
    for iFiles = 1:numel(filesPerControl)
        D = readfiles(filesPerControl(iFiles),param);
        %     Remove cells out of focus
        focus = getFocusFilteredData(D.data,param);
        D.data = D.data(focus,:);
        D.textdata = D.textdata(focus,:);
        allInten(cnt:cnt+size(D.data,1)-1,:) = D.data(:,intFeat);
        allD(cnt:cnt+size(D.data,1)-1,:) = D.data(:,param.datafeat);
        cnt = cnt+size(D.data,1);
    end
    if(cnt<mxRw)
        allD = allD(1:cnt-1,:);
        allInten = allInten(1:cnt-1,:);
    end
    ii = allInten>100;
    allD = allD(ii,:);
    allInten = allInten(ii,:);
    ii = allInten>quantile(allInten,.05) & allInten<quantile(allInten,.95);
    allD = allD(ii,:);
    allD = bsxfun(@minus,allD, param.meaninc(1,param.datafeat));
    allD = bsxfun(@rdivide,allD,sqrt(param.varinc(1,param.datafeat)));
    allD = allD(:,redFeatures);
%     allIndex = zeros(size(allD,1),1);
    newDensity= getDensityBasedSampling(allD,distanceType);
    allD = allD(newDensity,:);
    indx = phenograph( allD, 5);
    uIndx = unique(indx);
    for iIndx = 1:numel(uIndx)
        clsCentroidLevel1(bCnt,:) = mean(allD(indx==uIndx(iIndx),:));
        controlCategory(bCnt,1) = D.textdata(1,9);
        bCnt = bCnt+1;
    end
end
if(bCnt<(maxClsPerControl*numel(cntrl)))
    clsCentroidLevel1 = clsCentroidLevel1(1:bCnt-1,:);
    controlCategory = controlCategory(1:bCnt-1,:);
end
clear allD allInten D uIndx indx bCnt filesPerControl
%% Plot data

% scr = clsCentroidLevel1*cf(:,1:2);
% scr=compute_mapping(clsCentroidLevel1, 't-SNE', 2);
% plot(scr(:,1),scr(:,2),'or','MarkerFacecolor','r');
%% Plot the centroids obtained from above
numK = [2:2:60];
numUIndx = zeros(numel(numK),1);
for i = 1:numel(numK)
    indx = phenograph( clsCentroidLevel1, numK(i));
    numUIndx(i) = numel(unique(indx));
end
figure;plot(numK,numUIndx,'-*b');
%% Recluster centroids

indxL2 = phenograph( clsCentroidLevel1, 4);

%% Assign each data point to centroids
mxRws = 1000000;
allIndex = zeros(mxRws,1);
allTxt =cell(mxRws,1);
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
    ii = D.data(:,intFeat)>100;
    D.data = D.data(ii,param.datafeat);
    
    D.textdata = D.textdata(ii,9);
    D.data = bsxfun(@minus,D.data, param.meaninc(1,param.datafeat));
    D.data = bsxfun(@rdivide,D.data,sqrt(param.varinc(1,param.datafeat)));
    D.data = D.data(:,redFeatures);
    allIndex(cnt:cnt+size(D.data,1)-1,:) = knnsearch(clsCentroidLevel1,D.data);
    allTxt(cnt:cnt+size(D.data,1)-1,1) = D.textdata;
    cnt = cnt+size(D.data,1);
end
if(cnt<mxRws)
    allIndex = allIndex(1:cnt-1,:);    
    allTxt = allTxt(1:cnt-1,:);
end
nIdx = indxL2(allIndex);
clear D focus cnt mxRw

% Compute Level 2 centroids

uIndxL2 = unique(indxL2);
clsCentroidLevel2 = zeros(numel(uIndxL2),size(clsCentroidLevel1,2));
for i = 1:numel(uIndxL2)
    clsCentroidLevel2(i,:) = mean(clsCentroidLevel1(indxL2==uIndxL2(i),:));
end

% Compute distribution
uC = unique(allTxt);
uIndxL2 = unique(nIdx);
cDistribution = zeros(numel(uC),numel(uIndxL2));
for j = 1:numel(uC)
    jj = strcmpi(allTxt,uC{j,:});
    for i = 1:numel(uIndxL2)
        ii = nIdx==uIndxL2(i);        
        cDistribution(j,i) = sum(jj.*ii);
    end
end
D = pdist2(clsCentroidLevel2,clsCentroidLevel2,'euclidean');
[w,xst] = kruskal(1-eye(size(clsCentroidLevel2,1)),D);
% Compute t-SNE plot of centroids
scr=compute_mapping([clsCentroidLevel1;clsCentroidLevel2], 't-SNE', 2);
scr1 = scr(1:size(clsCentroidLevel1,1),:);
% Plot
figure;hold on;
for i = 1:numel(uC)
    ii = strcmpi(controlCategory,uC{i,:});
    plot(scr1(ii,1),scr1(ii,2),'o','MarkerFaceColor',param.maps(i,:),...
        'MarkerEdgeColor','none');
end
hold off;
legend(uC);title('Level 1 Centroids')


%%
D = pdist2(scr(size(clsCentroidLevel1,1)+1:end,:),scr(size(clsCentroidLevel1,1)+1:end,:),'euclidean');
[w,xst] = kruskal(1-eye(size(clsCentroidLevel2,1)),D);
h=viewMSTPie2(scr(size(clsCentroidLevel1,1)+1:end,:),...
    cDistribution',param.maps,uC,xst);
