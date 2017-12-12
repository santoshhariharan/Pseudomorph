function  multiBarHistOverlay( x,y,map,colName,ylabelVal )
%boxHistOverlay Creates a box plot with histogram overlay for every group
%in y
% Inputs:
% x: feature vector: 1 figure per column
% y: Groups
% Check number of inputs
if(nargin < 3)
    fprintf('@boxHistOverlay: Not enough inputs\n');
    return;
end
% Remove nan or Inf
rmoveRow=or(isnan(x),isinf(x));
x=x(~rmoveRow);
y=y(~rmoveRow);
% Check if y is a cell array of string
if(iscellstr(y))
    groupNames = unique(y);
    uy = unique(y);
    ny = uint8(zeros(numel(y),1));
    for i = 1:numel(uy)
        ny(strcmpi(uy{i,:},y),1) = i;
    end
    y = ny;
    uy = unique(y);
    clear ny;
else
    y = uint8(y);
    uy = unique(y);
    groupNames = cellstr(num2str([1:numel(uy)]'));
end
[m,n] = size(x);

if(m~=numel(y))
    fprintf('@boxHistOverlay: Dimension mismatch between x & y\n');
    return;
end

if(max(y)>size(map,1))
    map = jet(double(max(y)));
end
if(nargin==3)
    colName = cellstr(num2str(1:n));
    ylabelVal = 'Normalized Feature';
elseif(nargin==4)
    ylabelVal = 'Normalized Feature';
end

boxWidth = 1.5;
jitterVal = 0;
nBins = floor(sqrt(size(x,1)));
nBins = min([nBins 200]);

xStart = .5;
[~,xout] = hist(x(:,1),nBins);
ymin = quantile(x(:,1),.10);
ymax = quantile(x(:,1),.90);
randSamSize = 100000;
allLeg = [];
figure;hold on;
for iGroups = 1:numel(uy)
    ii = y==uy(iGroups);
%     x25 = quantile(x(ii,1),.25);
%     x75 = quantile(x(ii,1),.75);
%     xdata = [xStart xStart xStart+boxWidth xStart+boxWidth];
%     ydata = [x25 x75 x75 x25];
% %     patch(xdata,ydata,1,'Linewidth',.5,'Facecolor','none','Edgecolor',map(iGroups,:));
%     x1 = x(x>x75);x2 = x(x<x25);
%     if(numel(x1)>randSamSize)
%         x1 = x1(randperm(numel(x1),randSamSize));
%     end
%     if(numel(x2)>randSamSize)
%         x2 = x2(randperm(numel(x2),randSamSize));
%     end
%     plot(xStart+boxWidth/2+(jitterVal*randn(numel(x1),1)),x1,'o','Markersize',2,'MarkerFacecolor',map(iGroups,:),'MarkerEdgeColor','None');
%     plot(xStart+boxWidth/2+(jitterVal*randn(numel(x2),1)),x2,'o','Markersize',2,'MarkerFacecolor',map(iGroups,:),'MarkerEdgeColor','None');
    y1 = hist(x(ii,1),xout);
    y1 = y1./sum(y1);
    y1 = (y1-min(y1))./(max(y1)-min(y1));
%     y1 = y1*20;
    h1=barh(xout,y1+xStart,'ShowBaseLine','off','BaseValue',iGroups*5,'Edgecolor','none','Facecolor',map(iGroups,:));
%     line([xStart xStart],[ymin ymax],'Color',map(iGroups,:),'LineStyle','--','LineWidth',2);
    allLeg = [allLeg;h1];
    xStart = xStart+boxWidth+1;
end
hold off;
% ylim([-3 15]);
ylabel(ylabelVal,'Interpreter','None');
% set(gca,'XTickLabel','');
% legend(allLeg,groupNames,'Interpreter','None','Location','BestOutside');
% title(colName,'Interpreter','None');
end

