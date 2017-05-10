% Load sampled data
clear all;
clc;
pth='F:\Projects\Proteinlocalization\PseudoMorph\TestDataForPseudoMorph';
load(fullfile(pth,'sampleDataCellFeat.mat'));
load(fullfile(pth,'parameters.mat'));
% Perform feature reduction once with PCA based and once with unsupervised
% greedy optimization
[redFeatures, featureScores ]= unsupervisedPCASelect(data,30);
data=data(:,redFeatures);
% Community Detection
% Perform community detection over multiple values of K
% Use Jaccard distance & then euclidean distance
numNeighbours = [5:5:300];
numCls = zeros(numel(numNeighbours),1);
for jNeighbors = 1:numel(numNeighbours)
    indx = phenograph( data, numNeighbours(jNeighbors));
    numCls(jNeighbors) = numel(unique(indx));
end
clc;
% disp('Ready');
%%
setK = 50;
indx = phenograph( data, setK);
uIndx = unique(indx);
addpath(genpath('F:\Projects\Proteinlocalization\PseudoMorph\Code\drtoolbox'));
score = compute_mapping(data, 't-SNE', 2);
score = score(:,1:2);
rmpath(genpath('F:\Projects\Proteinlocalization\PseudoMorph\Code\drtoolbox'));
clusterCentroidPCA = zeros(numel(uIndx),size(score,2));
uC = unique(textdata);
cDistribution = zeros(numel(uC),numel(uIndx));
for j = 1:numel(uC)
    jj = strcmpi(textdata,uC{j,:});
    for i = 1:numel(uIndx)
        ii = indx==uIndx(i);
        clusterCentroidPCA(i,:) = mean(score(ii,:));
        cDistribution(j,i) = sum(jj.*ii);
    end
end
am = ones(numel(uIndx));
for i = 1:numel(uIndx)
    am(i,i) = 0;
end
D = pdist2(clusterCentroidPCA,clusterCentroidPCA,'euclidean');
[w,xst] = kruskal(am,D);

disp('Ready');
%% Plot results for viewing
h = figure;
plot(numNeighbours,numCls,'-o',...
        'MarkerEdgeColor','None','MarkerFaceColor','b');
    xlabel('#Neighbours(K)');ylabel('#Clusters');
    set(gca,'YTick',[5:5:75]);
savefig(h,fullfile(pth,'SampledDataClusterVNumNeighbors_All.fig'));

mp = jet(numel(uC));
h=viewMSTPie2(clusterCentroidPCA,cDistribution',mp,uC,xst);
savefig(h,fullfile(pth,'MSTPie_All.fig'));


h=figure; hold on;
for i = 1:numel(uC)
    ii=strcmpi(textdata,uC{i,:});
    plot(score(ii,1),score(ii,2),'o',...
        'MarkerEdgeColor','None','MarkerFaceColor',mp(i,:));
end
hold off;
xlabel('t-SNE1');ylabel('t-SNE2');zlabel('t-SNE3');
grid on;title('Sampled Data');
legend(uC);
savefig(h,fullfile(pth,'sampleDataFiguretSNE_All.fig'));

%% ROUND 2
% Perform community detection one control at a time
% Merge clusters by euclidean distance, by permutation test
% using wards method or KNN



