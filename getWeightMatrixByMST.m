function [ y,aPicked ] = getWeightMatrixByMST( mCent2,k,n )


if(nargin==1)
    n = 500;
    k = 10;
elseif(nargin==2)
    n = 500;
elseif(nargin<1)
    disp('@getWeightMatrixByMST: Incomplete arguments');
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
%     pp = pdist2(mCent2(idx,:),mCent2(idx,:));
    tmp = mCent2(idx,:);
    [kk,kkDis] = knnsearch(tmp,tmp,'k',k+1);
    kk = kk(:,2:end);
    kkDis = kkDis(:,2:end);
    pp = inf(size(tmp,1));
    for j = 1:size(kk,1)
        pp(j,kk(j,:)) =kkDis(j,:);
        pp(kk(j,:),j) =kkDis(j,:)';
    end
%     pp = exp(-pp);
%     pp(eye(size(pp,1))==1) = 0;
%     sumPP = sum(pp,2);
%     pp = bsxfun(@rdivide,pp,sumPP);
%     pp = 1-pp;
%     pp = (pp+pp');
    ad = 1-eye(numel(idx));
    [~,x,~] = kruskal(ad,pp);
    for j = 1:size(x,1)
        y(idx(x(j,1)),idx(x(j,2))) = y(idx(x(j,1)),idx(x(j,2)))+1;
        y(idx(x(j,2)),idx(x(j,1))) = y(idx(x(j,1)),idx(x(j,2)));
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
% ii = aPicked>=10;
% y = (y+y');
y =y./aPicked;
% aConnections(aConnections==0) = -inf;
% aConnections = 1-aConnections;


end

