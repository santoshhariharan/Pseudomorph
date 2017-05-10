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
cellAreaFeat = strcmpi('Ch1_MOR_Cell_area',param.datahdr);
filePrefix = '.txt';
% randPrc = .1;
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
fprintf('\nRemoved NAN\n');

% Remove incorrectly segmented & low intensity Objects
ii = allMorRatio <= .7;
allD = allD(ii,:);
allTxt = allTxt(ii,:);
allTxtOrg = allTxtOrg(ii,:);
allInten = allInten(ii,:);
% allMorRatio = allMorRatio(ii,:);
allMorIntensity = allMorIntensity(ii,:);
fprintf('#Cells Removed 4 intensity %i of %i\n',sum(~ii),numel(ii));

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
% Feature Selection/Reduction
% newHeader = param.datahdr(1,param.datafeat);
featureReduction = true;
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
minSamplePerControl = floor(.5*minSamplePerControl);
for i = 1:max(gps)
    ii = find(gps == i);    
    samIndex(ii(randperm(numel(ii),minSamplePerControl))) = true;
end
fprintf('minimum samples per control %i\n',minSamplePerControl);

% Cluster Samples data for high number of clusters

clc;
k = 5;
graphType = 'jaccard';
data4Clustering = allD(samIndex,:);
nGps = gps(samIndex);
mCent = nan(1000,size(allD,2));
mGrp = nan(1000,1);
cnt = 1;
for iControl = 1:max(gps)
    ii = nGps==iControl;
    fprintf('%s\n - %d\n',uControls{iControl,:},sum(ii));
    grpData = data4Clustering(ii,:);
    indx = phenograph(grpData,k,'graphtype',graphType);
    uIndx = unique(indx);
    for i = 1:numel(uIndx)
        mCent(cnt,:) = mean(grpData(indx==uIndx(i),:));
        mGrp(cnt,1) = iControl;
        cnt = cnt+1;
    end
end
disp('Done');
clear  grpData indx uIndx ii;
ii = sum(isnan(mCent),2)==0;
mCent = mCent(ii,:);
mGrp  = mGrp(ii,:);
m = size(mCent,1);
indx = knnclassify(data4Clustering, mCent, 1:m);
% indx = knnclassify(allD, mCent, 1:m);
uIndx = unique(indx);
nGps = gps(samIndex);
% gHatMat = false(numel(nGps),numel(uIndx));
% gMat = false(numel(nGps),max(nGps));
% for i = 1:max(nGps)
%     gMat(:,i) = nGps==i;
% end
% for i = 1:numel(uIndx)
%     gHatMat(:,i) = indx==uIndx(i);
% end
nij = zeros(numel(uIndx),max(nGps));
for i = 1:numel(uIndx)
    ii = indx==uIndx(i);
    for j = 1:max(nGps)
        jj = nGps ==j;
        nij(i,j) = sum(ii.*jj);
    end
end
% nij = double(gMat')*double(gHatMat);
% nij = nij';
ii = sum(nij,2)>0;
nij = nij(ii,:);
mCent = mCent(ii,:);
mGrp = mGrp(ii,:);
%
% nij = bsxfun(@rdivide,nij,sum(nij,1));
% Compute Enrichment
pp = enrich_score(nij);
enrichCutoff = .05./(numel(uIndx)*max(nGps));
ppLog = pp<enrichCutoff;
% Number of pure clusters
ppOrig = false(numel(uIndx),max(nGps));
for i = 1:max(nGps)
    iAll = false(1,max(nGps));
    iAll(1,i) = true;
    ii = sum(ppLog(:,~iAll),2) == 0;
    kk = sum(and(ppLog(:,iAll),ii));
    ppOrig(:,i) = and(ppLog(:,iAll),ii);
    fprintf('Group %i - %i\n',i,kk);
end

radMarker = zeros(numel(uIndx),1);
for i=1:numel(uIndx)
    radMarker(i) = sum(indx == uIndx(i));
end
radMarker = radMarker./sum(radMarker);
radMarker = ((radMarker - min(radMarker))./(max(radMarker)-min(radMarker)));
radMarker = floor(8*(radMarker) + 4);
clear pp ppLog i iAll ii kk;
clear gMat gHatMat
clear enrichCutoff k graphType data4Clustering
clear nGps allMorIntensity allMorRatio j jj i ii k
disp('Completed Mod 1');
%%
numDimsVis = 2;
dismat = bsxfun(@rdivide,nij,sum(nij,2));
pcaRedData = compute_mapping(mCent,'t-SNE',numDimsVis);
% pcaRedData = mdscale(pdist2(dismat,dismat),numDimsVis);
%%
% Find KNN graph
kVal = 20;
[idx, dis] = knnsearch(pcaRedData,pcaRedData,'K',kVal+1);
idx = idx(:,2:end);
% Plot data
figure; hold on;
% for iControl = 1:m
%     line([pcaRedData(iControl,1) pcaRedData(idx(iControl,:),1)'],...
%             [pcaRedData(iControl,2) pcaRedData(idx(iControl,:),2)'],...
%             'Color',[.5 .5 .5]);
% end

for iControl = 1:m
%     ii = mGrp == iControl;
    plot(pcaRedData(iControl,1),pcaRedData(iControl,2),'o','MarkerFaceColor',map(mGrp(iControl),:),...
        'MarkerEdgeColor','None','MarkerSize',radMarker(iControl));     
    
    
end

% set(gca,'XTick',[]);set(gca,'YTick',[]);
hold off; 
legend(uControls);

w = knn2jaccard(idx);
% w = idx2knn(idx,dis(:,2:end));
w = 1-full(w);
w = (w+w')/2;
scaleFact = 500;
adjMat = 1-eye(size(mCent,1));
viewScatterPie2(pcaRedData,indx,allTxt(samIndex,1),map,scaleFact);
[~,xst] = kruskal(adjMat, w);
viewMSTPie2(pcaRedData,nij,map,uControls,xst,scaleFact);

clear w adjMat a xst scaleFact idx dis kVal;
% w = pdist2(pcaRedData,pcaRedData,'cityblock');
% idx = idx(:,2:end);
% dis = dis(:,2:end);
% % [idx1, dis1] = knnsearch(dismat,dismat,'K',kVal+1);
% adjMat = false(size(mCent,1),size(mCent,1));
% w = zeros(size(mCent,1),size(mCent,1));
% d = zeros(size(mCent,1),size(mCent,1));
% for i = 1:size(idx,1)
%     adjMat(i,idx(i,:)) = true;
%     w(i,idx(i,:)) = dis(i,:);
% %     d(i,idx1(i,:)) = dis1(i,:);
% end
% w = computePathlength( adjMat, w );

% w = pdist2(mCent,mCent,'cityblock');
% d = pdist2(dismat,dismat,'cityblock');

% adjMat = adjMat+adjMat';
% adjMat = adjMat>0;
% w = bsxfun(@rdivide,w,sum(w,2));
% d = bsxfun(@rdivide,d,sum(d,2));
% w = w+d;w = w+w';
%
% figure; hold on;
% for i = 1:numel(uIndx)
%     if(numDimsVis==3)
%         plot3(pcaRedData(i,1),pcaRedData(i,2),pcaRedData(i,3),'o',...
%             'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','None',...
%             'MarkerSize',radMarker(i));
%         line([pcaRedData(i,1) pcaRedData(idx(i,:),1)'],...
%             [pcaRedData(i,2) pcaRedData(idx(i,:),2)'],...
%             [pcaRedData(i,3) pcaRedData(idx(i,:),3)'],'Color',[.5 .5 .5]);
%     else
%         plot(pcaRedData(i,1),pcaRedData(i,2),'o',...
%             'MarkerFaceColor',[0 0 1],'MarkerEdgeColor','None',...
%             'MarkerSize',3);
%         line([pcaRedData(i,1) pcaRedData(idx(i,:),1)'],...
%             [pcaRedData(i,2) pcaRedData(idx(i,:),2)'],...
%             'Color',[.5 .5 .5]);
%     end
% end
% hold off;