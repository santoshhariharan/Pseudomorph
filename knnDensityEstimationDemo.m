
% Read & sample 20% of the data
% This shoudl be accomplished for every control separately
% Keep an index for the total number of data points
% For now load everything, later only write out x% of the data points
% initially

clear all;close all;clc;
numelGP2 = 50;
data = [2*randn(1000,2); (7+3.*randn(numelGP2,2))];
grps = ones(size(data,1),1);
grps(end-numelGP2:end) = 2;
% clear fid userSelectedControls;
%% Plot original data
figure;hold on;
plot(data(grps==1,1),data(grps==1,2),'*r');
plot(data(grps==2,1),data(grps==2,2),'*b');
hold off;title('Simulated gaussian data');
%% Take random samples
p = randperm(size(data,1),.2*size(data,1));
cellstoKeep = false(numel(grps),1);
cellstoKeep(p) = true;
pData = data(cellstoKeep,:);
pGrps = grps(cellstoKeep,:);
figure;hold on;
plot(pData(pGrps==1,1),pData(pGrps==1,2),'*r');
plot(pData(pGrps==2,1),pData(pGrps==2,2),'*b');
hold off;title('Randomly Sampled data');


%% Get density estimate for every point
allLdensity = zeros(numel(cellstoKeep),1);
% distancePref = .5;% Median of distances
% numNeighbours =  120;
% [~,D] = knnsearch(pData,pData,'K',numNeighbours,'Distance','Euclidean');
% idx = idx(:,2:end);
% D = D(:,2:end);
%%
distancePref = .5;% Median of distances
x =pdist2(pData,pData);
distanceThreshold = quantile(x(:),distancePref);
% distanceThreshold = distanceThreshold.*(numNeighbours./size(pData,1));
% dThresh = D<=distanceThreshold;
lDensity = sum(x<=distanceThreshold,2)./size(pData,1);
lDensity(lDensity==1)=.99;
lDensity(lDensity==0)=.01;
allLdensity(cellstoKeep) = lDensity;
% figure;histogram(lDensity);xlabel('Local density');title('Training LD histogram');
figure;scatter(pData(:,1),pData(:,2),30,lDensity,'filled');colorbar;
title('Density of random samples')
%% Load the rest of the data to get estimate of local density
% userSelectedControls = getappdata(0,'alldata');
% fid = fopen(filename,'r');
% data = fread(fid,[size(userSelectedControls.textdata,1) sum(param.datafeat)],'double');
% fclose(fid);
%%
% data = data(~cellstoKeep,:);
numNeighbours = 10;
[idx,D] = knnsearch(data(cellstoKeep,:),data(~cellstoKeep,:),'K',...
                    numNeighbours,'Distance','Euclidean');
%%
% D = 1./D;
D = exp(-D);
pp = bsxfun(@rdivide,D,sum(D,2));
ff = reshape(lDensity(idx,:),size(idx,1),size(idx,2));
cc = sum(pp.*ff,2);
allLdensity(~cellstoKeep) = sum(pp.*ff,2);

% figure;histogram(allLdensity);xlabel('Local density');title('All data LD histogram');
disp('@@@');
figure;scatter(data(:,1),data(:,2),30,allLdensity,'filled');colorbar;
title('Local Density Plot');
% return;
%%
newallLdensity = 1-allLdensity;
newDensity = newallLdensity>rand(numel(allLdensity),1);
% figure;histogram(newallLdensity(newDensity));
%%
figure;scatter(data(newDensity,1),data(newDensity,2),30,allLdensity(newDensity),'filled');
colorbar;
title('Sampled points based on LD');
%%
dataForClustering = data(newDensity,:);
[idx] = knnsearch(dataForClustering,dataForClustering,'K',50,'Distance','Euclidean');
idx = idx(:,2:end);
sim = knn2jaccard(idx);
[com] = cluster_jl(sim,1,0);
[~,maxMod] = max(com.MOD);
indx = com.COM{1,maxMod}';
uIndx = unique(indx);
s=['*r';'*b';'*g';'*k'];
figure;hold on;
for i = 1:numel(uIndx)
    plot(dataForClustering(indx==uIndx(i),1),dataForClustering(indx==uIndx(i),2),s(i,:));
end
hold off;
return;
%%
if(sum(newDensity)>10000)
    newDensity = newDensity
    nBin = [0:.05:1];
    for i = 1:numel(nBin)
        
    end
end
data2 = data(newDensity,:);
[idx] = knnsearch(data2,data2,'K',numNeighbours,'Distance','Euclidean');
idx = idx(:,2:end);
sim = knn2jaccard(idx);
[com] = cluster_jl(sim,1,0);
[~,maxMod] = max(com.MOD);
indx = com.COM{1,maxMod}';    
numClusters = numel(unique(indx));
disp('####');














