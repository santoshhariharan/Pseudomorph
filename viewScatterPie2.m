function [ h ] = viewScatterPie2(x,y,yhat,map,scalingFact)
%viewScatterPie Plots data in x using grouping in y with colours of yhat
% 
% 

if(size(x,2)>2)
    error('@viewScatterPie: Only two dmensional plot allowed');    
end

if(size(y,1)~=size(yhat,1))
    error('@viewScatterPie: y and yhat must be of same length');
end

uYhat = unique(yhat);
if(size(map,1) ~= numel(uYhat))
    map = jet(numel(uYhat)); % Default Color map
end
uY = unique(y);
% Compute distribution matrix
disMat = zeros(numel(uYhat),numel(uY));
for iTreatments = 1:numel(uYhat)
    ii = strcmpi(uYhat{iTreatments,:},yhat);
    for jClusters = 1:numel(uY)
        jj = y==uY(jClusters);
        disMat(iTreatments,jClusters) = sum(ii.*jj);
    end
end

rsize = sum(disMat,1);
rsize = rsize./sum(rsize);

% scalingFact = 60;
% rsize = rsize*scalingFact;
% First plot clusters with lines from centroid of each cluster
rangeX = max(x) - min(x);
maxRSize = min(rangeX/(2*scalingFact));
if(maxRSize<0)
    maxRSize = 2;
end
rsize = (rsize - min(rsize))./(max(rsize) - min(rsize));
rsize = (maxRSize-1)*rsize +1;

    
h=figure;hold all;
% for iClusters = 1:numel(uY)
%     ii = y == uY(iClusters);
%     if(sum(ii)==0)
%         continue;
%     end
%     plot(x(ii,1),x(ii,2),'o',...
%         'MarkerSize',4,'MarkerEdgeColor','none',...
%         'MarkerFaceColor',[.6 .6 .6]);
%     cCenter = repmat(x(uY(iClusters),:),sum(ii),1);
%     x1 = [cCenter(:,1) x(ii,1)];
%     y1 = [cCenter(:,2) x(ii,2)];
%     line(x1',y1','Color',[.6 .6 .6],'MarkerSize',6,'MarkerEdgeColor','none',...
%         'MarkerFaceColor',[.6 .6 .6]);
% end
% hold off;
%
% In the same figure plot data using patches - Uses inbuilt Matlab function
% logic
% 
% numSmallPatch = 30;
numTreatments = size(disMat,1);
% angleStart = 0;
numSteps = 30;
% rsize = 5;
% hold on;
for iClusters = 1:numel(uY)
    startAngle = 0;  
    treatProportion = disMat(:,iClusters)./sum(disMat(:,iClusters));
%     treatProportion = bsxfun(@rdivide,disMat(:,iClusters),sum(disMat,2));
%     treatProportion = treatProportion./sum(treatProportion);
    cCenter = x(uY(iClusters),:);
    hall=[];
    for jTreatment =  1:numTreatments
        step = treatProportion(jTreatment)*2*pi/numSteps;
        theta = [0 0:step:treatProportion(jTreatment)*2*pi 0];
        theta = theta + startAngle;
        [x1, y1] = pol2cart(theta,[0 rsize(iClusters).*(ones(1,numel(theta)-2)) 0]);
        x1 = x1 + cCenter(1,1);
        y1 = y1 + cCenter(1,2);
        h1 = patch(x1,y1,jTreatment,'FaceColor',map(jTreatment,:),'EdgeColor','none');
        hall = [hall;h1];
        startAngle = max(theta);  
       
    end    
%     legend(uYhat)
end

hold off;
legend(hall,uYhat);
% clear hall

end

