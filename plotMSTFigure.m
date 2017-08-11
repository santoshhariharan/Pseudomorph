function  plotMSTFigure( pcaRedData,xst,mGrp,map )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
m = size(xst,1);
figure; hold on;
for i = 1:m
    line([pcaRedData(xst(i,1),1) pcaRedData(xst(i,2),1)'],...
            [pcaRedData(xst(i,1),2) pcaRedData(xst(i,2),2)'],...
            'Color',[.7 .7 .7]);
end
uMgrp = unique(mGrp);
for iControl = 1:numel(uMgrp)
%     if(sum(rmGrp == iControl)>0)
%         continue;
%     end
    
    ii = mGrp == uMgrp(iControl);
    plot(pcaRedData(ii,1),pcaRedData(ii,2),'o','MarkerFaceColor',map(iControl,:),...
        'MarkerEdgeColor','None','MarkerSize',3);      
end
% legend(uControls);
% set(gca,'XTick',[]);set(gca,'YTick',[]);
hold off;

end

