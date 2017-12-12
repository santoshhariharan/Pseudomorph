function [ h ] = viewMSTPie2(x,disMat,map,xst,lineWidth,scalingFact,radius)
%viewMSTPie2 Plots data in x using grouping in y with colours of yhat
% Plots cluster centroids based on x. 
% Input:
% x - cluster centroids #clusters x 2
% map - map for treatment distribution of clusters
% disMat - #Treatments x#clusters
% names - Names of treatment
% xst - Minimum Spanning tree

if(nargin==4)
    lineWidth = .1*ones(size(xst,1),1);
    scalingFact = .1;
    radius = ones(size(x,1),1);
end
if(size(x,2)>2)
    error('@viewScatterPie: Only two dmensional plot allowed');    
end

if(size(map,1)<size(disMat,2))
    error('@viewScatterPie: Map does not include all treatments');
end

axisRange = abs(max(x)-min(x));
lineWidth = lineWidth*5;
% Radius of the highest is 10% of axis range
maxRadius = max(5*scalingFact*axisRange);
minRadius = min(.5*scalingFact*axisRange);

% rsize = sum(disMat,2);
radius = radius./sum(radius);
% radius = (radius - min(radius))./(max(radius)-min(radius));
rsize = abs(maxRadius - minRadius)*radius + minRadius;

% scalingFact = max(rsize);

% scalingFact = 50;
% scalingFact = rMaxWidth./scalingFact;
% rsize = rsize*scalingFact;

% Resize scaling between 4 & 20
% maxR = 10;
% minR = 4;
% rsize = ((maxR-minR)*((rsize - min(rsize))/(max(rsize)-min(rsize)))+minR);
% First plot clusters with lines from centroid of each cluster
% rsize = floor(rsize);


% In the same figure plot data using patches - Uses inbuilt Matlab function
% logic
% 
% numSmallPatch = 30;
numTreatments = size(disMat,2);
% angleStart = 0;
numSteps = 30;
% rsize = 5;
% hold on;
h=figure;hold all;
for i = 1:size(xst,1)
    coord = x(xst(i,:),:)';    
    line(coord(1,:),coord(2,:),'Linewidth',lineWidth(i),'Color',[.7 .7 .7]);    
end

for iClusters = 1:size(disMat,1)
    startAngle = 0;  
    treatProportion = disMat(iClusters,:)./sum(disMat(iClusters,:));
    cCenter = x(iClusters,:);
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



hold off;axis tight
% legend(hall,names);
% clear hall

end

