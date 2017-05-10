function sampledPoints = getDensityBasedSampling(data,distanceType,oulierdensity)
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
x = x(:,1:5);
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
pointDensity = 1-pointDensity;
sampledPoints = pointDensity>rand(numel(pointDensity),1);
sampledPoints(pointDensity<=quantile(pointDensity,oulierdensity)) = false;
fprintf('Initial: %.0f Final %.0f\n',numel(pointDensity),sum(sampledPoints));

% cellstokeep = false(size(data,1),1);
% % Take 30% of data @ random
% cellstokeep(randperm(size(data,1),round(.3.*size(data,1)))) = true;
% x = pdist2(data(cellstokeep,:),data(cellstokeep,:),distanceType);
% x = sort(x,2,'ascend');
% x = x(:,2:end);
% medDis = mean(x(:));
% lDensity = (sum(x<medDis,2)-1)./(sum(cellstokeep)-1);
% allLdensity(cellstokeep,1) = lDensity;
% [idx,D] = knnsearch(data(cellstokeep,:),data(~cellstokeep,:),'K',...
%     numneighbors,'Distance',distanceType);
% % D = exp(-D);
% allLdensity(~cellstokeep) = sum(bsxfun(@rdivide,D,sum(D,2)).*reshape(lDensity(idx,:),...
%     size(idx,1),size(idx,2)),2);
% clear D;
% 
% allLdensity = 1-allLdensity;
% sampledPoints = allLdensity>rand(numel(allLdensity),1);
end