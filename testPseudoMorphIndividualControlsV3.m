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

clear;clc;
fprintf('Starting pseudomorph\n');
pth='F:\Projects\Proteinlocalization\PseudoMorph\TestDataForPseudoMorph';
load(fullfile(pth,'parameters.mat'));% Load parameter file
param.rootpath = pth;
% cntrl = {'CB5-ER';'ER-PRO';'ERGIC';'FL-VAMP5';'FL-VAMP2';...
%     'GOLGI-GT';'GOLGIN';'METAXIN';...
%     'PQC-PSS1';'RAB5A';'RAB7A';'MAO';'MITO-CCO';...
%     'TABIK'};
intensityFeature = 'Ch2_INT_Cell_intensity';
intFeat = strcmpi(intensityFeature,param.datahdr);
filePrefix = 'Controls';
randPrc = .1;
numFeatures = 30;
fNames = dir(pth);
columnForControls = 9;
featureReduction = true;
clustersPerLandmark = true;
%% Module 2 - Compute mean and variance
% fprintf('Module 2.......\n');
% md = computeMeanVarianceCls(pth,param);
% param.meaninc = md.meaninc;
% param.varinc = md.varinc;
% param.maxColVal = md.maxColVal;
% param.minColVal = md.minColVal;
% clear md;clc;
% fprintf('Completed Mean & Variance Computation\n');


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
allInten = allInten(jj,:);
param.meaninc = mean(allD);
param.varinc = var(allD);
allD = zscore(allD);
allInten = zscore(allInten);
clear D focus cnt tok iFiles mxRw allInten
%%
featureReduction = false;
if(featureReduction)
    numFeatures = 50;
    redFeatures = unsupervisedGreedyFS(allD,numFeatures);
else
    redFeatures = true(1,size(allD,2));
end
% [cf]= princomp(allD(:,redFeatures));
allD = allD(:,redFeatures);



%% Find number of samples per control & cap it to 10000

uControls = unique(allTxt);
numSamplesPerControl = 10000; % For clustering
for i = 1:numel(uControls)
    numSamplesPerControl = min(numSamplesPerControl,sum(strcmpi(uControls{i,:},allTxt)));
end


%% Module 4: Per control
fprintf('Module 4.......\n');
valueOfK = 20;
if(clustersPerLandmark)
    maxClsPerControl = 30;
    distanceType = 'Euclidean';
    uControls = unique(allTxt);
    bCnt = 1;
    clsCentroidLevel1 = zeros(maxClsPerControl.*numel(uControls),sum(redFeatures));
    controlCategory = cell(maxClsPerControl.*numel(uControls),1);
    elementsPerCluster  = zeros(maxClsPerControl.*numel(uControls),1);
    for iCnt = 1:numel(uControls)
        
        ii = find(strcmpi(uControls{iCnt,:},allTxt));
        %     numSamplesPerControl = numel(ii);
        jj = randperm(numel(ii),numSamplesPerControl);
        
        %     allIndex = zeros(size(allD,1),1);
        %     newDensity= getDensityBasedSampling(allD,distanceType);
        %     allD = allD(newDensity,:);
        dataForClustering = allD(ii(jj),:);
        indx = phenograph( dataForClustering, valueOfK,'distance','Euclidean',...
            'graphtype','euclidean');
        uIndx = unique(indx);
        fprintf('%s\t%i\n',uControls{iCnt,:},numel(uIndx));
        for iIndx = 1:numel(uIndx)
            fprintf('   C%i\t%i\n',iIndx,sum(indx==uIndx(iIndx)));
            clsCentroidLevel1(bCnt,:) = median(dataForClustering(indx==uIndx(iIndx),:));
            elementsPerCluster(bCnt) = sum(indx==uIndx(iIndx));
            controlCategory(bCnt,1) = uControls(iCnt,:);
            bCnt = bCnt+1;
        end
        
    end
    if(bCnt<(maxClsPerControl*numel(uControls)))
        clsCentroidLevel1 = clsCentroidLevel1(1:bCnt-1,:);
        controlCategory = controlCategory(1:bCnt-1,:);
        elementsPerCluster = elementsPerCluster(1:bCnt-1,:);
    end
else
    sampleData = nan(numSamplesPerControl.*numel(uControls),size(allD,2));
    sampleTxt = cell(numSamplesPerControl.*numel(uControls),1);
    cnt = 1;
    for iCnt = 1:numel(uControls)
        ii = find(strcmpi(uControls{iCnt,:},allTxt));
        %     numSamplesPerControl = numel(ii);
        jj = randperm(numel(ii),numSamplesPerControl);
        sampleData(cnt:cnt+numel(jj)-1,:) = allD(ii(jj),:);
        sampleTxt(cnt:cnt+numel(jj)-1,:) = allTxt(ii(jj),:);
        cnt = cnt+numel(jj);
    end
    ii = sum(isnan(sampleData),2)==0;
    sampleData = sampleData(ii,:);
    sampleTxt = sampleTxt(ii,:);
    indx = phenograph( sampleData, valueOfK,'distance','Euclidean',...
            'graphtype','euclidean');
    
    uIndx = unique(indx);
    fprintf('%s\t%i\n',uControls{iCnt,:},numel(uIndx));
    clsCentroidLevel1 = nan(numel(uIndx),size(allD,2));
    elementsPerCluster = nan(numel(uIndx),1);
    for iIndx = 1:numel(uIndx)
        fprintf('   C%i\t%i\n',iIndx,sum(indx==uIndx(iIndx)));
        clsCentroidLevel1(iIndx,:) = median(sampleData(indx==uIndx(iIndx),:));
        elementsPerCluster(iIndx) = sum(indx==uIndx(iIndx));        
    end
    ii = sum(isnan(clsCentroidLevel1),2)==0;
    clsCentroidLevel1 = clsCentroidLevel1(ii,:);
    elementsPerCluster = elementsPerCluster(ii,1);
end
clear D uIndx indx bCnt dataForClustering ii jj
%% Plot data

% scr = clsCentroidLevel1*cf(:,1:2);
% scr=compute_mapping(clsCentroidLevel1, 't-SNE', 2);
% plot(scr(:,1),scr(:,2),'or','MarkerFacecolor','r');
%% Plot the centroids obtained from above
% numK = [2:2:60];
% numUIndx = zeros(numel(numK),1);
% for i = 1:numel(numK)
%     indx = phenograph( clsCentroidLevel1, numK(i));
%     numUIndx(i) = numel(unique(indx));
% end
% figure;plot(numK,numUIndx,'-*b');


% scr1 = scr(1:size(clsCentroidLevel1,1),:);
%% Plot Filter with less than 50
ii = elementsPerCluster >=50;
mSize = 100*ones(sum(ii),1);


        
clsCentroidLevel1 = clsCentroidLevel1(ii,:);
controlCategory = controlCategory(ii,:);
elementsPerCluster = elementsPerCluster(ii,:);
tmp = ((.8)*((elementsPerCluster - min(elementsPerCluster))./...
            (max(elementsPerCluster) - min(elementsPerCluster))) + .1);
mSize = ((.8)*((elementsPerCluster - min(elementsPerCluster))./...
            (max(elementsPerCluster) - min(elementsPerCluster))) + .1).*mSize;
        
scr=compute_mapping(clsCentroidLevel1, 't-SNE', 3);
figure;hold on;
for i = 1:numel(uControls)
    ii = strcmpi(controlCategory,uControls{i,:});
    %     plot3(scr(ii,1),scr(ii,2),scr(ii,3),'o','MarkerFaceColor',param.maps(i,:),...
    %         'MarkerEdgeColor','none');
    scatter3(scr(ii,1),scr(ii,2),scr(ii,3),mSize(ii,1),param.maps(i,:),...
                    'filled');
end
hold off;
legend(uControls);title('Level 1 Centroids');
return;

%%
%% Recluster centroids

indxL2 = phenograph( clsCentroidLevel1, 2,'distance','Euclidean',...
                    'graphtype','euclidean');
% Compute Level 2 centroids

uIndxL2 = unique(indxL2);
clsCentroidLevel2 = zeros(numel(uIndxL2),size(clsCentroidLevel1,2));
for i = 1:numel(uIndxL2)
    clsCentroidLevel2(i,:) = median(clsCentroidLevel1(indxL2==uIndxL2(i),:));
end
scr=compute_mapping(clsCentroidLevel2, 'PCA', 3);
figure;
plot3(scr(:,1),scr(:,2),scr(:,3),'o','MarkerFaceColor',param.maps(i,:),...
        'MarkerEdgeColor','none');

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
