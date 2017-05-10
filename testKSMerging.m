load('parameters.mat');
filename = fullfile('F:\Projects\Proteinlocalization\PseudoMorph\TestDataForPseudoMorph',...
                '_controls_temp.bin');
treatcol = 9;
gpnames = {'GOLGIN'};
userSelectedControls = getappdata(0,'alldata'); % Get Metadata information 
fid = fopen(filename,'r');
data = fread(fid,[size(userSelectedControls.textdata,1) sum(param.datafeat)],'double');
fclose(fid);
% Get randomsamples based on controls
grps = getGroupIndices(userSelectedControls.textdata(:,treatcol),...
                                gpnames);
cellstoKeep = false(numel(grps),1);
for i =1:numel(unique(grps))
    f = find(grps == i);
    p = randperm(numel(f),10000);
    cellstoKeep(f(p),1) = true;
end
data = data(cellstoKeep,:);
grps = (grps(cellstoKeep,1));
redFeatures = unsupervisedGreedyFS(data,80);
data = data(:,redFeatures);
%%
p = randperm(sum(cellstoKeep),floor(sum(cellstoKeep)/2));
grps(p,1) = 2;
clusterCount= 1;
clusterIndex = zeros(sum(cellstoKeep),1);
for iGroups = 1:max(grps)    
    idx = knnsearch(data(grps==iGroups,:),data(grps==iGroups,:),'K',125,'Distance','Euclidean');
    idx = idx(:,2:end);
    sim = knn2jaccard(idx);
    [com] = cluster_jl(sim,1,0);
    [~,maxMod] = max(com.MOD);
    indx = com.COM{1,maxMod}';
    indx = indx+clusterCount+1;
    clusterCount = max(indx)+clusterCount;
    clusterIndex(grps==iGroups) = indx;
    fprintf('#Group: %i #Clusters: %d\n',iGroups,numel(unique(indx)));
end
%%
% clc;
k = -1:-1:-10;
% alphaVal = [.05 1e-2];
numC = zeros(numel(k),1);
for i = 1:numel(k)
    fprintf('#Iter: %i\n',i);
    newIdx = mergeClustersKSTest( data,clusterIndex, 10.^k(i), 'Euclidean' );
    numC(i) = numel(unique(newIdx)); 
%     re
end
fprintf('#Cls-Grp1: %i #Cls-Grp2: %i\n',numel(unique(clusterIndex(grps==1))),numel(unique(clusterIndex(grps==2))));
% figure;hold on;
% semilogx(10.^k,repmat(numel(unique(clusterIndex(grps==1))),numel(k),1),'-xb');
% semilogx(10.^k,repmat(numel(unique(clusterIndex(grps==2))),numel(k),1),'-xb');
% semilogx(10.^k,numC,'-xr');







