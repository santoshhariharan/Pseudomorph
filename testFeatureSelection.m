% Compare feature selection types
% Compare random forests, SDA, unsupervised greedy method

% load paramater files
compareClassification = true;
load('parameters.mat');
filename = fullfile('F:\Projects\Proteinlocalization\PseudoMorph\TestDataForPseudoMorph',...
                '_controls_temp.bin');
treatcol = 9;
gpnames = {'ERGIC';'GOLGIN';'TABIK'};
% read data and metadata            
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

if(~compareClassification) % Create new data set for clustering
%     Get Means
    cellstoKeep = false(size(data,1),1);
    for i =1:numel(unique(grps))
        f = find(grps==i);
        mD = median(data(f,:));
        p = knnsearch(data(f,:),mD,'K',100);
        cellstoKeep(f(p),1) = true;
%         disp('@')
    end
    data = data(cellstoKeep,:);
    grps = (grps(cellstoKeep,1));
end
[~,sc] = princomp(data);
figure;hold on;
plot(sc(grps==1,1),sc(grps==1,2),'*r');
plot(sc(grps==2,1),sc(grps==2,2),'*b');
plot(sc(grps==3,1),sc(grps==3,2),'*g');
hold off;
clear fid cellstoKeep userSelectedControls filename
%% Test feature selection on random data
clc;
k = [5:5:size(data,2)];
% y = unsupervisedGreedyFS(data,20);
classAccuMean = zeros(numel(k),1);
randAccuMean = zeros(numel(k),1);
classAccuStd = zeros(numel(k),1);
randAccuStd = zeros(numel(k),1);
for i = 1:numel(k)
    
%     selectedFeatures = false(1,size(data,2));
    redFeatures = unsupervisedGreedyFS(data,k(i));
    
%     selectedFeatures(1,redFeatures) = true;
    if(compareClassification)
        [ac,rc] = runRFclassifier(data,grps,redFeatures);
    else
        [ac,rc] = runClustering(data,grps,redFeatures);
    end
    classAccuMean(i) = mean(ac);
    randAccuMean(i) = mean(rc);
    classAccuStd(i) = std(ac);
    randAccuStd(i) = std(rc);
    fprintf('#Features: %i  #Accu: %.2f  #RandAccu: %.2f\n',...
                    sum(redFeatures),classAccuMean(i),randAccuMean(i));
end

figure;hold on;
errorbar(k,classAccuMean*100,classAccuStd*100,'-xr');
errorbar(k,randAccuMean*100,randAccuStd*100,'-ob');
hold off;
xlabel('Number of features');ylabel('Classification Accuracy');
title('Feature selection comparison');
legend({'Greedy Unsupervised';'Random'});





