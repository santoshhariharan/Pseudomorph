function [ indx, mIndx, objectAssign ] = consensusClustering( x,options,plotting)
%consensusClustering Perform consunsus clustering by creating ensembles
%                   from data
% Ref: 'Cluster Ensembles – A Knowledge Reuse Framework for Combining 
%           Multiple Partitions'
% Input:
% x - Data (m observations, n variables)
% options - A structure with the following fields
%   .distance - Distance to be computed - default 'Jaccard'
%   .algorithm - clsutering algorithm: choice between, 'AP', phenograph,
%                   k-means
%   .pref - Preference for AP - default 'pmin6'
%   .kcls - Value of k for kmeans
%   .k - Value for NN graph for phenograph
%   .sampling - 'random' or 'densitybased' - default 'random'
%   .repeats - Number of ensembles - default '50'
%   .grouped - Grouping variable - m x1 vector indicatin the groups. If
%           more than 1 group, clustering is done separately for each group and
%           combined later
% Output:
% indx - Final index of clsutering
% mIndx - Index clusering for every ensemble
% 
if(nargin < 1 || isempty(x))
    error('Not enough inputs');
end

[m,n] = size(x);


% Check the number of arguments

if(nargin == 1)
    options.distance = 'Jaccard';
    options.algorithm = 'phenograph';    
    options.pref = 'pmin6';
    options.kcls = 10;
    options.k = 70;
    options.sampling = 'random';
    options.repeat = 10;
    options.grouped = ones(m,1);
    plotting = false;
    
elseif(nargin == 2)
    plotting = false;
end

% Get minimum sample number
minSamplePerGroup = 1000;
uGrp = unique(options.grouped);
for i = 1:numel(uGrp)
    minSamplePerGroup = min(minSamplePerGroup,sum(options.grouped == uGrp(i)));
end
fprintf('#Groups: %.0f\t minsample %.0f \n',numel(uGrp),minSamplePerGroup);
% Start groupwise clustering
% indx = uint8(zeros(m,1));
mIndx = uint8(zeros(m,options.repeat));
hGraph = false(m,500);clsCount = 1;
allCentroids = [];
for iRpt = 1:options.repeat % Start number of repeats 
    mCentroids = nan(100,n);cnt = 1;
    for iGrp = 1:1
        ii = options.grouped == uGrp(i);
        dataPerGroup = x(ii,:);
        if(strcmpi(options.sampling,'densitybased'))
            jj = true(sum(ii),1);
        else
            jj = randperm(sum(ii),minSamplePerGroup);
        end
        dataPerGroup = dataPerGroup(jj,:);
        groupIndex = phenograph( dataPerGroup, options.k,'distance',options.distance); % Perform clsutering
        uGIndx = unique(groupIndex);
        for iCentroids = 1:numel(uGIndx)
            mCentroids(cnt,:) = mean(dataPerGroup(groupIndex == uGIndx(iCentroids),:));
            cnt = cnt+1;            
        end        
    end
    jj = sum(isnan(mCentroids),2) == 0;
    mCentroids = mCentroids(jj,:);
    allCentroids = [allCentroids;mCentroids];
%     Assign clusters by closest euclidean distance    
    mIndx(:,iRpt) = knnclassify(x,mCentroids,uint8([1:size(mCentroids,1)]'),1);
%     Create Hypergraph
    for iCls = 1:size(mCentroids,1)
        hGraph(:,clsCount) = mIndx(:,iRpt) == iCls;
        clsCount = clsCount+1;
    end
end
clear dataPerGroup
fprintf('Starting hypergraph clustering\n');
if(options.repeat == 1)
    indx = mIndx;
    objectAssign = mIndx;    
    return;
end

% For plotting purposes
[~,allCentroids] = princomp(allCentroids);
% allCentroids = compute_mapping(allCentroids,'t-SNE',3);
% Combine different clusterings 
% Create Hypergraph
hGraph = hGraph(:,sum(hGraph,1) >0);

% Recluster the hypergraph
% Compute jaccard similarity

hGraph = hGraph';
[p,~] = size(hGraph);
[idx, dis] = knnsearch(hGraph,hGraph,'K',options.k+1,...
                'Distance','Jaccard');
idx = idx(:,2:end);
dis = dis(:,2:end);
sim = nan(p*(p-1)/2,3);cnt = 1;
for i = 1:size(hGraph,1)
    sim(cnt:cnt+options.k-1,1) = i;
    sim(cnt:cnt+options.k-1,2) = idx(i,:)';
    sim(cnt:cnt+options.k-1,3) = 1-dis(i,:)';
%     for j = i+1:size(hGraph,1)        
%         sim(cnt,1)  = j;sim(cnt,2) = i;
%         cnt = cnt +1;
%     end
    cnt = cnt+options.k;
end
ii = sum(isnan(sim),2)==0;
sim = sim(ii,:);
% s = 1-pdist(hGraph,'Jaccard');
% ii = s~=0;
% s = s(:,ii)';
% sim = sim(ii',:);
s = sim(:,3);
% Plot nearest neighbour Graph in PCA
if(plotting)
    b = size(sim,1);
    figure;hold on;
    plot3(allCentroids(:,1),allCentroids(:,2),allCentroids(:,3),'o','Markerfacecolor',[.6 .6 .6],...
        'MarkerEdgecolor','none','Markersize',4);
    for i = 1:b
        tmp1 = allCentroids(sim(i,1),:);
        tmp2 = allCentroids(sim(i,2),:);
        line([tmp1(1,1) tmp2(1,1)],[tmp1(1,2) tmp2(1,2)],[tmp1(1,3) tmp2(1,3)],'Color',[.4 .4 .4]);
    end
    hold off;
    title('Nearest neighbour graph');
end
G = sparse(sim(:,1),sim(:,2),s,p,p);
u = triu(G,1);
v = tril(G,-1);
G = v + u';
labels = phenograph( x, 1, 'G',G );
uLabels = unique(labels);
objectAssign = zeros(numel(uLabels),m);
if(plotting)
    mp = jet(numel(uLabels));
    figure;hold on;
    for i = 1:numel(uLabels)
        ii = labels == uLabels(i);
        objectAssign(i,:) = mean(hGraph(ii,:));
        plot3(allCentroids(ii,1),allCentroids(ii,2),allCentroids(ii,3),'o','Markerfacecolor',mp(i,:),...
            'MarkerEdgecolor','none','Markersize',4);
    end
    hold off;
    title('Combined labels');
end
% objectAssign = objectAssign';
[~,indx] = max(objectAssign,[],1);

objectAssign = objectAssign';

% Plot centroids

fprintf('Completed Consesus clustering');




end

