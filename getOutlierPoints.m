function oulierPoints = getOutlierPoints(data,distanceType,oulierdensity)
% getDensityBasedSampling: Computes the local density of each data point
% and samples based on inverse of the local density. If the local density
% is less than outlier density, those points are ignored
% Inputs:
% data - data martaix m x n
% distancetype - Limit by pdist
% oulierdensity - Minimum density of points 
% numneighbors - Number of neighbors for knn search
% sampledPoints = zeros(size(data,1),1);
if(nargin == 2)
    oulierdensity = .01;
end

I = randperm(size(data,1),round(.1*size(data,1)));
x = pdist2(data(I,:),data(I,:),distanceType);
x = sort(x,2,'ascend');
x = x(:,2:end);
% x = x(:,1:5);
medDis = quantile(x(:),.5);

% Compute Distance between data points
pointDensity = zeros(size(data,1),1);
fprintf('Computing local Density.........')
for i = 1:size(data,1)
    if(pointDensity(i) == 0)
%         d = pdist2(data(i,:),data,distanceType);
        d = bsxfun(@minus,data,data(i,:));
        d = sqrt(sum(d.*d,2));
        pointDensity(i) = sum(d<medDis);
        pointDensity(d<.75*medDis) = pointDensity(i); % Smoothing
    end    
    if mod(i,1000)==1 || i == size(data,2)
        fprintf('\b\b\b\b\b%4d%%',floor(i/length(pointDensity)*100));
    end
end
fprintf('\n');
pointDensity = pointDensity./length(pointDensity);
oulierPoints = pointDensity<oulierdensity;
fprintf('#Percent Outlier Points %f\n',sum(oulierPoints)*100/numel(oulierPoints));
end