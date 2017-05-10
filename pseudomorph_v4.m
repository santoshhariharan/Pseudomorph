% Load paramateres
clc;clear;
warning off;
pth = 'F:\Projects\Proteinlocalization\PseudoMorph\TestDataForPseudoMorph';
load(fullfile(pth,'parameters.mat'));
dimReduction = true;
if(~isfield(param,'meaninc'))
    md = computeMeanVarianceCls(pth,param);
    param.meaninc = md.meaninc;
    param.varinc = md.varinc;
    param.maxColVal = md.maxColVal;
    param.minColVal = md.minColVal;
    clear md;clc;
    fprintf('Completed Mean & Variance Computation');
    % Save parameters
    save(fullfile(pth,'parameters.mat'),'param','-append');
else
    disp('Skipping mean variance computation');
end

if(exist('ClusterOutput','file')==7)
    rmdir(fullfile(pwd,'ClusterOutput'),'s');
end
% Map controls to filenames
% cntrl = {'ergic';'golgin';'tabik';'ER-PRO';'Golgi-GT'};
cntrl = {'golgin';'ergic';'tabik'};
% cntrl = {'ergic'};
mapFilesToControls;

%% perform feature reduction
if(dimReduction)
    unsupervisedFeatureReduction;
else
    finalFeatures = true(1,sum(param.datafeat));
end
disp('Completed Unsupervised selection');
%%
dimCellsThreshold = 100;% Minimum value of intensity for cell to be considred for analysis
% intensityLowerThreshold = .01;
% intensityUpperThreshold = .99;
% channelName = 'Ch2_INT_Cell_intensity'; % Intensity channel name;
intensityChannel = strcmpi(param.datahdr,'Ch2_INT_Cell_intensity');
distanceType='Euclidean';
% Go each control read correspoding files and perform community detection
clusterFileMap =  cell(1000,3);
mkdir('ClusterOutput');
clusterCount = 0;
numNeighbours = [20];
maxSamples = 150000;
%%
numClustersNeighbors = zeros(numel(numNeighbours),1);
numClusterPerControl = zeros(numel(numNeighbours),size(cntrl,1));
for jNeighbors = 1:numel(numNeighbours)
    clusterCount=0;maxSamples = 150000;
for iControls = 1:size(cntrl,1)
    fprintf('Control name: %s\n',cntrl{iControls,:});
    files = fileMap(strcmpi(fileMap(:,2),cntrl{iControls,:}),1);
    D = readfiles(files,param);
%     allLdensity = zeros(size(D.data,1),1);
    intensityLowerThreshold = quantile(D.data(:,intensityChannel),.01);
    intensityUpperThreshold = quantile(D.data(:,intensityChannel),.99);
    
    D.textdata = D.textdata(D.data(:,intensityChannel)>dimCellsThreshold,:);
    D.data = D.data(D.data(:,intensityChannel)>dimCellsThreshold,:);
    iii = and(D.data(:,intensityChannel)>intensityLowerThreshold,...
        D.data(:,intensityChannel)<intensityUpperThreshold);
    D.textdata = D.textdata(iii,:);
    D.data = D.data(iii,:);
%     D.textdata = D.textdata(D.data(:,intensityChannel)<intensityUpperThreshold,:);
%     D.data = D.data(D.data(:,intensityChannel)<intensityUpperThreshold,:);
    focus = getFocusFilteredData(D.data,param);
    D.textdata = D.textdata(focus,:);
    D.data = D.data(focus,:);
    D.data = bsxfun(@minus,D.data, param.meaninc);
    D.data = bsxfun(@rdivide,D.data,sqrt(param.varinc));
    data = D.data(:,param.datafeat);
    data = data(:,finalFeatures);
    textdata = D.textdata;
    clear D;
%     if(size(data,1)>maxSamples) 
%         indToKeep  = 
%         data = data(randperm(size(data,1),maxSamples),:);
%     end
    newDensity= getDensityBasedSampling(data,distanceType);
%     [idx] = knnsearch(data(newDensity,:),data(newDensity,:),...
%         'K',numNeighbours,'Distance','Euclidean');
%     idx = idx(:,2:end);
%     sim = knn2jaccard(idx);
    indx = phenograph( data(newDensity,:), numNeighbours(jNeighbors), 'Distance','Euclidean' );
%     [com] = cluster_jl(sim,1,0);
%         com = cluster_jl_cpp(sim,1,0,0,1);
%     [~,maxMod] = max(com.MOD);
%     indx = com.COM{1,maxMod}';
    numClusters = numel(unique(indx));
    numClusterPerControl(jNeighbors,iControls) = numClusters;
    fprintf('--#Clusters: %i\n',numClusters);
    allIndex = zeros(size(data,1),1);
    allIndex(newDensity) = indx;
    %     Compute centroids and assign the rest of the data to centroids
    clusterCenter = zeros(numClusters,size(data,2));
    ii = find(newDensity);
    uIdx = unique(indx);
    clustersize = zeros(numel(uIdx),1);
    for i = 1:numClusters
        clusterCenter(i,:) = mean(data(ii(indx == uIdx(i)),:));
        clustersize(i) = numel(ii(indx == uIdx(i)));
    end
    
    % Assign points to nearest cluster
    allIndex(~newDensity) = knnclassify(data(~newDensity,:),...
        clusterCenter, uIdx,1);
    
%     Remove clusters with very low number of cells
%     clusterCenter = clusterCenter(clustersize>.001*size(data,1));
%     clustersize = clustersize(clustersize>.001*size(data,1));
%   Write data to file per cluster
    numClusters = size(clusterCenter,1);
    for i = 1:numClusters
        filename = ['Cluster_' cntrl{iControls} '_' num2str(i+clusterCount) '.bin'];
        txtFileName = ['Cluster_' cntrl{iControls} '_' num2str(i+clusterCount) '_text.txt'];
        filename = fullfile(pwd,'ClusterOutput',filename);
        txtFileName= fullfile(pwd,'ClusterOutput',txtFileName);
        clustersize(i) = sum(allIndex==uIdx(i));
        try
        fid = fopen(filename,'w');       
        fwrite(fid,data(allIndex==uIdx(i),:),'double');
        fclose(fid);
        writestr(txtFileName,textdata(allIndex==uIdx(i),:),'overwrite');
        catch e
            fclose(fid);
            rethrow(e);
        end
        clusterFileMap{i+clusterCount,1} = filename;
        clusterFileMap{i+clusterCount,2} = clustersize(i);
        clusterFileMap{i+clusterCount,3} = cntrl{iControls,:}; 
    end
    clusterCount = clusterCount+numClusters;
end
clear clusterCount numClusters 
clear allIndex allLdensity indx clusterCenter clusterCount
clear clustersize com data  textdata
clear fid filename files i iControls idx 
clear ii indx maxMod
% clear distanceType numNeighbours
clear newDensity numFeatures  sim tmpFeat uIdx
save(fullfile(pwd,'ClusterOutput','clusterParameters.mat'),'clusterFileMap','finalFeatures');
mergeClusters;
numClustersNeighbors(jNeighbors) = numel(unique(mergingIndex(:,iIterations+1)));
end
figure;plot(numNeighbours,numClustersNeighbors,'o','Markerfacecolor','b');
clear distanceType numNeighbours
% return;
%% Merge Clusters

mergeClusters;
viewdata;