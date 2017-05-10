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
    p = randperm(numel(f),5000);
    cellstoKeep(f(p),1) = true;
end

data = data(cellstoKeep,:);
grps = (grps(cellstoKeep,1));
redFeatures = unsupervisedGreedyFS(data,80);
data = data(:,redFeatures);
clear fid cellstoKeep userSelectedControls filename;
%%
k = [5:5:145];
numClusters = zeros(numel(k),1);
idx = knnsearch(data,data,'K',k(end)+1,'Distance','Euclidean');
idx = idx(:,2:end);
for i = 1:numel(k)
    sim = knn2jaccard(idx(:,1:k(i)));
    [com] = cluster_jl(sim,1,0);
    [~,maxMod] = max(com.MOD);
    indx = com.COM{1,maxMod}';    
    numClusters(i) = numel(unique(indx));
    fprintf('#K: %i #Clusters: %d\n',k(i),numClusters(i));
end
figure;
plot(k,numClusters,'-xb');
xlabel('Number of neighbours');ylabel('Number of clusters');


