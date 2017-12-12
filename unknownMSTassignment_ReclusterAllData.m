clear; clc;
% Load the centroid file
%%
centroidFileName = 'centroidPerControl_Feat30F_5K.mat';
unKnowmRootDir = 'F:\Projects\Proteinlocalization\PseudoMorph\Unknowns';

load(centroidFileName);
reclusterInitialCentroids = true;
if(reclusterInitialCentroids)
    k = 5;
    indx = phenograph(mCent,k,'graphtype','Jaccard');
    uIndx = unique(indx);
    mCent2 = zeros(numel(uIndx),size(mCent,2));
    mFrac2 = zeros(numel(uIndx),1);
    for j = 1:numel(uIndx)
        jj = find(indx == uIndx(j));
        mCent2(j,:) = mean(mCent(jj,:));
        mFrac2(j,1) = numel(jj)/numel(indx);        
    end
else
    mCent2 = mCent;
    mGrp2 = mGrp;
    mFrac2 = mFraction;
end
clc;

%% Compute Ensemble MST
mstCutoffValue= .5;
aConnections = getWeightMatrixByMST(mCent2);
% aConnections = tril(aConnections,-1);
jj = aConnections<=mstCutoffValue;
aConnections(jj) = 0;

%% Load unknowns with filters

fprintf('Starting pseudomorph\n');
% pth2paramfile='F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data';
load(fullfile(unKnowmRootDir,'parameters.mat'));% Load parameter file
param.rootpath = unKnowmRootDir;
intensityFeature = 'Ch2_INT_Cell_intensity';
intFeat = strcmpi(intensityFeature,param.datahdr);
nucAreaFeat = strcmpi('Ch1_MOR_Nucleus_area',param.datahdr);
cellAreaFeat = strcmpi('Ch1_MOR_Cytoplasm_area',param.datahdr);
filePrefix = '.txt';
fNames = dir(unKnowmRootDir);
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

% Remove Loe intensity Objects
ii =allInten>=100; % Value of this needs optimization
allD = allD(ii,:);
allTxt = allTxt(ii,:);
allTxtOrg = allTxtOrg(ii,:);
allInten = allInten(ii,:);
allMorRatio = allMorRatio(ii,:);
% allMorIntensity = allMorIntensity(ii,:);

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
% load('redFeatures_60.mat');
allD = allD(:,redFeatures);
clear D focus cnt tok iFiles mxRw 
clear allMorRatio
clear ii jj kk rho mxRw 
clear intensityFeature intFeat roiIntFeat nucAreaFeat cellAreaFeat
clear fNames filePrefix randPrc

allD = bsxfun(@minus,allD,meanD);
allD = bsxfun(@rdivide,allD,stdD);
%% % Assign unknown to MST

[row, column] = find(aConnections>0);
dis = disPoint2liine(allD,mCent2(row,:),mCent2(column,:));
I = find(aConnections>0);
% wght = aConnections(:);
% wght = wght(I);
% dis = 1./dis;
% wght = repmat(wght',size(dis,1),1);
% dis = dis.*wght;
[nDis,I] = sort(dis,2,'ascend');
x = [row column];
% I = I(:,1);

%% Plot tSNE plots for each


redData = compute_mapping(mCent2,'t-SNE',2);
figure; hold on;
for i = 1:max(mGrp)
    ii = mGrp2 == i;
%     plot(redData(ii,1),redData(ii,2),'o',...
%         'MarkerFaceColor',map(i,:),'MarkerEdgeColor','none',...
%         'MarkerSize',floor(mFrac2(ii,1)*20));
    scatter(redData(ii,1),redData(ii,2),floor(mFrac2(ii,1)*500),...
        map(i,:),'filled'); 
%     scatter3(redData(ii,1),redData(ii,2),redData(ii,3),floor(mFrac2(ii,1)*500),...
%         map(i,:),'filled');
end
hold off;  
% axis off;
%

figure; hold on;
for i = 1:size(x,1)
    line([redData(x(i,1),1) redData(x(i,2),1)],[redData(x(i,1),2) redData(x(i,2),2)],...
        'color',[.7 .7 .7],'linewidth',aConnections(x(i,1),x(i,2))*3);
end
for i = 1:max(mGrp)
    ii = mGrp2 == i;
%     plot(redData(ii,1),redData(ii,2),'o',...
%         'MarkerFaceColor',map(i,:),'MarkerEdgeColor','none',...
%         'MarkerSize',floor(mFrac2(ii,1)*20));
    scatter(redData(ii,1),redData(ii,2),floor(mFrac2(ii,1)*500),...
        map(i,:),'filled'); 
end
hold off;
% axis off;

%% Project points onto Line
y = x(I(:,1),:);
[~,pdis] = projectPointonLine( allD,mCent2(y(:,1),:),mCent2(y(:,2),:) );
pdis(pdis>1)=1;
% Compute New Coordinates
newCoord = zeros(size(allD,1),size(redData,2));
for i = 1:size(allD,1)
    C1 = redData(y(i,1),:);
    C2 = redData(y(i,2),:);
    C3 = C2-C1;
    cSquared = sqrt(dot(C3,C3));
    newCoord(i,:) = C1 + (pdis(i).*C3);
end
%% Plot on Graph

figure; hold on;
for i = 1:size(x,1)
    line([redData(x(i,1),1) redData(x(i,2),1)],[redData(x(i,1),2) redData(x(i,2),2)],...
        'color',[.7 .7 .7],'linewidth',aConnections(x(i,1),x(i,2))*3);
end
for i = 1:max(mGrp)
    ii = mGrp2 == i;
%     plot(redData(ii,1),redData(ii,2),'o',...
%         'MarkerFaceColor',map(i,:),'MarkerEdgeColor','none',...
%         'MarkerSize',floor(mFrac2(ii,1)*20));
    scatter(redData(ii,1),redData(ii,2),floor(mFrac2(ii,1)*1000),...
        map(i,:),'filled'); 
end

% for i = 1:size(newCord,1)
    scatter(newCoord(:,1),newCoord(:,2),20,...
        [1 1 1],'filled');
% end
hold off;axis off;
return;
%% Per Group
m = numel(mGrp2);
grp2Retain = [1 2 3 6 7 11];
idx2Keep = false(m,1);
for i = 1:numel(grp2Retain)
    idx2Keep = or(idx2Keep,mGrp2 == grp2Retain(i));
end
mCentTmp = mCent2(idx2Keep,:);
aConnections = getWeightMatrixByMST(mCentTmp);
mGrpTmp = mGrp2(idx2Keep);
nRedData = redData(idx2Keep,:);
mFracTmp = mFrac2(idx2Keep);
[r,c] = find(aConnections>.5);
x1 = [r c];
figure; hold on;
for i = 1:size(x1,1)
    line([nRedData(x1(i,1),1) nRedData(x1(i,2),1)],[nRedData(x1(i,1),2) nRedData(x1(i,2),2)],...
        'color',[.7 .7 .7],'linewidth',aConnections(x1(i,1),x1(i,2))*5);
end
for i = 1:max(mGrpTmp)
%     if(sum(grp2Retain==i)==0)
%         continue;
%     end
    ii = mGrpTmp == i;
%     plot(redData(ii,1),redData(ii,2),'o',...
%         'MarkerFaceColor',map(i,:),'MarkerEdgeColor','none',...
%         'MarkerSize',floor(mFrac2(ii,1)*20));
    scatter(nRedData(ii,1),nRedData(ii,2),floor(mFracTmp(ii,1)*1000),...
        map(i,:),'filled'); 
end
hold off;axis off;


