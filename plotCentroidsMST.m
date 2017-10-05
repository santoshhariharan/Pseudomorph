file2load = 'F:\Projects\Proteinlocalization\PseudoMorph\Code\centroidPerControl_RedFeat_5K.mat';

redData = compute_mapping(mCent,'t-SNE',2);

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
