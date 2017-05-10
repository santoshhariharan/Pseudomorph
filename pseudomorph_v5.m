
%% Load paramateres
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
% cntrl = {'golgin'};
% cntrl = {'ergic'};
fileMap = mapFilesToControls(pth,cntrl);
manualFiltering = true;
modeVal = false; % If mode is true, the clustering is done per control else clustering
             % is done on all data
testMode = true;
%% Choose Features for removing noisy data
% param=rmfield(param,'intensityLowerThreshold');
intensityChannel = strcmpi(param.datahdr,'Ch2_INT_Cell_intensity');
morphChannel = or(strcmpi(param.datahdr,'Ch1_MOR_Cell_area'),...
    strcmpi(param.datahdr,'Ch1_MOR_Cytoplasm_area'));
if(~isfield(param,'intensityLowerThreshold'))
    %         Read data from text file
%     if(modeVal)% Per control analysis
        intensityLowerThreshold = zeros(numel(cntrl),1);
        intensityUpperThreshold = zeros(numel(cntrl),1);
        morphometricFilterLowerThreshold = zeros(numel(cntrl),1);
        morphometricFilterUpperThreshold = zeros(numel(cntrl),1);
        for iControls = 1:size(cntrl,1)
            files = fileMap(strcmpi(fileMap(:,2),cntrl{iControls,:}),1); 
            [data,~] = readfilesMod(files,param);
            data = data(data(:,intensityChannel)>100,:);
            intensityLowerThreshold(iControls) = quantile(data(:,intensityChannel),.05);
            intensityUpperThreshold(iControls) = quantile(data(:,intensityChannel),.95);            
            d = data(:,morphChannel);
            d = d(:,2)./d(:,1);
            morphometricFilterLowerThreshold(iControls) = quantile(d,.05);
            morphometricFilterUpperThreshold(iControls) = quantile(d,.95);
%             ii = (d>quantile(d,.05) & d<quantile(d,.95));
        end
        param.intensityLowerThreshold = median(intensityLowerThreshold);
        param.intensityUpperThreshold = median(intensityUpperThreshold);
        param.morphometricFilterLowerThreshold = median(morphometricFilterLowerThreshold);
        param.morphometricFilterUpperThreshold = median(morphometricFilterUpperThreshold);
        param.intensityChannel = intensityChannel;
        param.morphChannel = morphChannel;
%         param.dimCellsThreshold = 100;
        save(fullfile(pth,'parameters.mat'),'param','-append');
%     end
fprintf('Complated Lower & Upper intnesity threshold computation\n')
end
% Print
if(param.intensityLowerThreshold<100)
    param.intensityLowerThreshold = 100;
end

fprintf('Lower: %4.3f\nUpper: %4.3f\n',param.intensityLowerThreshold,...
                param.intensityUpperThreshold);

%% perform feature reduction/selection
if(dimReduction)
    finalFeatures = unsupervisedFeatureReduction(fileMap,cntrl,param,'pcaselect');
else
    finalFeatures = true(1,sum(param.datafeat));
end
disp('Completed Unsupervised selection');
%%
dimCellsThreshold = 100;% Minimum value of intensity for cell to be considred for analysis
% intensityLowerThreshold = .01;
% intensityUpperThreshold = .99;
% channelName = 'Ch2_INT_Cell_intensity'; % Intensity channel name;


distanceType='Euclidean';
% Go each control read correspoding files and perform community detection
clusterFileMap =  cell(1000,3);
mkdir('ClusterOutput');
clusterCount = 0;
numNeighbours = [20:20:100];
% maxSamples = 150000;
%%
numClustersNeighbors = zeros(numel(numNeighbours),1);
numClusterPerControl = zeros(numel(numNeighbours),size(cntrl,1));
for jNeighbors = 1:numel(numNeighbours)
    clusterCount=0;maxSamples = 150000;
    for iControls = 1:size(cntrl,1)
        fprintf('Control name: %s\n',cntrl{iControls,:});
        if(modeVal)
            files = fileMap(strcmpi(fileMap(:,2),cntrl{iControls,:}),1);
        else
            files = fileMap(:,1);
            iControls = size(cntrl,1)+1;
        end
        files = files(randperm(size(files,1)),:);
        if(jNeighbors ==1)
            [data,textdata] = readfilesMod(files,param);
            %     Z score data - Normalization
            data = bsxfun(@minus,data, param.meaninc);
            data = bsxfun(@rdivide,data,sqrt(param.varinc));
            data = data(:,param.datafeat);
            data = data(:,finalFeatures);
            %     Get density based sampling
            newDensity= getDensityBasedSampling(data,distanceType);
        end
        
        indx = phenograph( data(newDensity,:), numNeighbours(jNeighbors));
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
        c = computeDistribution(allIndex,textdata(:,9));
        continue;
        %     Remove clusters with very low number of cells
        %     clusterCenter = clusterCenter(clustersize>.001*size(data,1));
        %     clustersize = clustersize(clustersize>.001*size(data,1));
        %   Write data to file per cluster
        numClusters = size(clusterCenter,1);
        
        for i = 1:numClusters
            if(modeVal)
                filename = ['Cluster_' cntrl{iControls} '_' num2str(i+clusterCount) '.bin'];
                txtFileName = ['Cluster_' cntrl{iControls} '_' num2str(i+clusterCount) '_text.txt'];
            else
                filename = ['Cluster_' num2str(i+clusterCount) '.bin'];
                txtFileName = ['Cluster_' num2str(i+clusterCount) '_text.txt'];
            end
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
        if(~modeVal)
            %         break;
        end
    end
    clear clusterCount
    clear allIndex indx clusterCenter clusterCount
    clear clustersize com
    clear fid filename files i iControls idx
    clear ii indx maxMod
    % clear distanceType numNeighbours
    clear  numFeatures  sim tmpFeat uIdx
    
    save(fullfile(pwd,'ClusterOutput','clusterParameters.mat'),'clusterFileMap','finalFeatures');
    if(modeVal)
        mergeClusters;
        numClustersNeighbors(jNeighbors) = numel(unique(mergingIndex(:,iIterations+1)));
    else
        numClustersNeighbors(jNeighbors) = numClusters;
    end
end
figure;plot(numNeighbours,numClustersNeighbors,'o','Markerfacecolor','b');
clear distanceType numNeighbours
return;
%% Merge Clusters

% mergeClusters;
viewdata;