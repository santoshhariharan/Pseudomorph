
% Read & sample 10% of the data
% This shoudl be accomplished for every control separately
% Keep an index for the total number of data points
% For now load everything, later only write out x% of the data points
% initially
clc;clear all;
featureReduction = true;
% reductionRatio= .1;
load('parameters.mat');
filename = fullfile('F:\Projects\Proteinlocalization\PseudoMorph\TestDataForPseudoMorph',...
                '_controls_temp.bin');
treatcol = 9;
gpnames = {'golgin'};
userSelectedControls = getappdata(0,'alldata'); % Get Metadata information 
try
fid = fopen(filename,'r');
data = fread(fid,[size(userSelectedControls.textdata,1) sum(param.datafeat)],'double');
fclose(fid);
catch exception
    fclose(fid);
    rethrow(exception)
end
numSamples = floor(.1*size(data,1));
fprintf('#UniformSamples: %d\n',numSamples);
% numSamples = 5000;
% Get randomsamples based on controls
grps = getGroupIndices(userSelectedControls.textdata(:,treatcol),...
                                gpnames);
cellstoKeep = false(numel(grps),1);
for i =1:numel(unique(grps))
    f = find(grps == i);
    p = randperm(numel(f),numSamples);
    cellstoKeep(f(p),1) = true;
end
if(featureReduction)
    redFeatures = unsupervisedGreedyFS(data(cellstoKeep,:),80);
else
    redFeatures = true(size(data,2));
end
clear fid userSelectedControls;
% Get density estimate for every point
allLdensity = zeros(numel(cellstoKeep),1);
x = pdist2(data(cellstoKeep,redFeatures),data(cellstoKeep,redFeatures));
%%
distancePref = .5;% Median of distances
distanceThreshold = quantile(x(x>0),distancePref);
dThresh = x<distanceThreshold;
lDensity = sum(x<distanceThreshold,2)./numSamples;
lDensity(lDensity==1)=.99;
lDensity(lDensity==0)=.01;
allLdensity(cellstoKeep) = lDensity;
figure;histogram(lDensity);xlabel('Local density');title('Training LD histogram');
% [~,score] = princomp(data(cellstoKeep,redFeatures));
% score = score(:,1:3);
% figure;scatter3(score(:,1),score(:,2),score(:,3),30,lDensity,'filled');colorbar;
% xlabel('PC1');ylabel('PC2');zlabel('PC3');
% title('PCA plot showing density');
clear x;
%% Load the rest of the data to get estimate of local density
% userSelectedControls = getappdata(0,'alldata');
[idx,D] = knnsearch(data(cellstoKeep,redFeatures),data(~cellstoKeep,redFeatures),'K',...
                    10,'Distance','Euclidean');
D = exp(-D);
allLdensity(~cellstoKeep) = sum(bsxfun(@rdivide,D,sum(D,2)).*reshape(lDensity(idx,:),...
                                size(idx,1),size(idx,2)),2);

figure;histogram(allLdensity);xlabel('Local density');title('All data LD histogram');
disp('@@@');
%% Sample data based on density
allLdensity = 1-allLdensity;
newDensity = allLdensity>rand(numel(allLdensity),1);
figure;histogram(allLdensity(newDensity));
xlabel('Local density');title('Sampled data based on inverse LD histogram');
%% Perfrom community detection
numNeighbours = 200;
fprintf('#NumNeighbors: %d\n',numNeighbours);
fprintf('#NumSamplesAfterDensitySampling: %d\n',sum(newDensity));
allIndex = zeros(size(data,1),1);
[idx] = knnsearch(data(newDensity,redFeatures),data(newDensity,redFeatures),...
                'K',numNeighbours,'Distance','Euclidean');
idx = idx(:,2:end);
sim = knn2jaccard(idx);
% [com] = cluster_jl(sim,1,0);
com = cluster_jl_cpp(sim,1,0);
[~,maxMod] = max(com.MOD);
indx = com.COM{1,maxMod}';    
numClusters = numel(unique(indx));
allIndex(newDensity) = indx;
disp('####');
%
% allIndex(~newDensity) = knnclassify(data(~newDensity,redFeatures),...
%                         data(newDensity,redFeatures), indx,numNeighbours);
% Find cluster centers
clusterCenter = zeros(numClusters,sum(redFeatures));
ii = find(newDensity);
uIdx = unique(indx);
for i = 1:numClusters    
    clusterCenter(i,:) = mean(data(ii(indx == uIdx(i)),redFeatures));
end

% Assign points to nearest cluster
allIndex(~newDensity) = knnclassify(data(~newDensity,redFeatures),...
                        clusterCenter, uIdx,1);
%% Plot data
uIndx = unique(allIndex);
numElementsPerCluster = zeros(numel(uIndx),1);
for i = 1:numel(uIndx)
    numElementsPerCluster(i) = sum(allIndex==uIndx(i));
    fprintf('#Cluster: %i #Elements %d\n',i,numElementsPerCluster(i));
end

figure;
semilogy(1:numel(uIndx),numElementsPerCluster,'or','MarkerFaceColor','r');
hold on;
semilogy(1:numel(uIndx),repmat(100,1,numel(uIndx)),'-b');
semilogy(1:numel(uIndx),repmat(400,1,numel(uIndx)),'-k');hold off;
ylabel('Number of cells')
disp('Completed Assignment');
