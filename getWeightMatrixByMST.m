function [ y ] = getWeightMatrixByMST( mCent2,n )


if(nargin==1)
    n = 500;
end
m = size(mCent2,1);
y = zeros(m);
aPicked = zeros(m);
numPointsPerTry = floor(.7*m);
fprintf('MST........................');
for i = 1:n
    fprintf('\b\b\b%3.0f',i);
    idx = randperm(m,numPointsPerTry);
    aPicked(idx,idx) = aPicked(idx,idx)+1;
    pp = pdist2(mCent2(idx,:),mCent2(idx,:));
    ad = 1-eye(numel(idx));
    [~,x,~] = kruskal(ad,pp);
    for j = 1:size(x,1)
        y(idx(x(j,1)),idx(x(j,2))) = y(idx(x(j,1)),idx(x(j,2)))+1;
    end
%     figure;hold on;
%     plot(mCent2(:,1),mCent2(:,2),'o','MarkerFacecolor',...
%         [.2 .2 .2],'MarkerEdgeColor','none');
%     
%     for k = 1:size(x,1)
%     line([mCent2(idx(x(k,1)),1) mCent2(idx(x(k,2)),1)],[mCent2(idx(x(k,1)),2) mCent2(idx(x(k,2)),2)],...
%         'color',[.5 .5 .5],'linewidth',2,'Marker','o','MarkerFacecolor','g',...
%         'MarkerEdgeColor','none');
%     end
%     hold off;
%     axis off;
end
y =y./aPicked;
% aConnections(aConnections==0) = -inf;
% aConnections = 1-aConnections;
y = y+y';

end

