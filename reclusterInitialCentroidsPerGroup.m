
clear; 
%%
clc;
reclusterData = true;
file2load = 'F:\Projects\Proteinlocalization\PseudoMorph\Code\centroidPerControl_RedFeat_5K.mat';
load(file2load);
if(reclusterData)
    k = 20;
    mCent2 = zeros(50,size(mCent,2));
    mGrp2 = zeros(50,1);
    mFrac2  = zeros(50,1);
    cnt = 1;
    clsCnt = 0;
    mIndx = zeros(size(mCent,1),1);
    for i = 1:max(mGrp)
        ii = find(mGrp == i);
        indx = phenograph(mCent(ii,:),k,'graphtype','Jaccard');
        indx = indx + clsCnt;
        clsCnt = clsCnt+max(indx)+1;
        mIndx(ii,1) = indx;
        uIndx = unique(indx);
        for j = 1:numel(uIndx)
            jj = find(indx == uIndx(j));
            mCent2(cnt,:) = mean(mCent(ii(jj),:));
            mGrp2(cnt,1) = i;
            mFrac2(cnt,1) = numel(jj)/numel(indx);
            cnt = cnt+1;
        end
    end
    cnt  =cnt-1;
    if(cnt<50)
        mCent2 = mCent2(1:cnt,:);
        mGrp2 = mGrp2(1:cnt,1);
        mFrac2 = mFrac2(1:cnt,1);
    end
else
    mCent2 = mCent;
    mGrp2 = mGrp;
    mFrac2 = mFraction;
end
clc;
% fprintf(%i\n)
for i = 1:max(mGrp2)
    fprintf('%i\n',sum(mGrp2==i));
end
disp('#####');
%% Plot Sammon plots for each
redData = compute_mapping(mCent2,'t-SNE',2);
%%
figure; hold on;
for i = 1:max(mGrp2)
    ii = mGrp2 == i;
%     plot(redData(ii,1),redData(ii,2),'o',...
%         'MarkerFaceColor',map(i,:),'MarkerEdgeColor','none',...
%         'MarkerSize',floor(mFrac2(ii,1)*20));
    scatter(redData(ii,1),redData(ii,2),floor(mFrac2(ii,1)*500),...
        map(i,:),'filled'); 
end
hold off;  
axis off;
%% Get Ensemble MST
clc;
aConnections = getWeightMatrixByMST(mCent2);
[r,c] = find(aConnections>0);
% Compute Euclidean Distance & MST
% pp = pdist2(mCent2,mCent2);
ad = 1-eye(size(mCent2,1));
% [~,x,~] = kruskal(ad,aConnections);
x = [r c];
figure; hold on;
for i = 1:size(x,1)
    line([redData(x(i,1),1) redData(x(i,2),1)],[redData(x(i,1),2) redData(x(i,2),2)],...
        'color',[.7 .7 .7],'linewidth',aConnections(x(i,1),x(i,2))*5);
end
for i = 1:max(mGrp2)
    ii = mGrp2 == i;
%     plot(redData(ii,1),redData(ii,2),'o',...
%         'MarkerFaceColor',map(i,:),'MarkerEdgeColor','none',...
%         'MarkerSize',floor(mFrac2(ii,1)*20));
    scatter(redData(ii,1),redData(ii,2),floor(mFrac2(ii,1)*500),...
        map(i,:),'filled'); 
end
hold off;axis off;
%% Plot different Pathways
% ER - ERGIC - Golgi - Lysosomes
m = numel(mGrp2);
grp2Retain = [  4 5 6 7 8];
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
    scatter(nRedData(ii,1),nRedData(ii,2),floor(mFracTmp(ii,1)*500),...
        map(i,:),'filled'); 
end
hold off;axis off;
%% Put Text
textVal = cellstr(num2str([1:numel(mGrp2)]'));
figure; hold on;
for i = 1:max(mGrp2)
    ii = mGrp2 == i;
%     plot(redData(ii,1),redData(ii,2),'o',...
%         'MarkerFaceColor',map(i,:),'MarkerEdgeColor','none',...
%         'MarkerSize',floor(mFrac2(ii,1)*20));
    scatter(redData(ii,1),redData(ii,2),floor(mFrac2(ii,1)*500),...
        map(i,:),'filled'); 
    text(redData(ii,1),redData(ii,2),textVal(ii,1),'Fontsize',10,'Color','k');
    
end
hold off;axis off;
