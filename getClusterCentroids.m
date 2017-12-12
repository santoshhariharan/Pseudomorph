function [ xCent,f ] = getClusterCentroids( x,idx,statType )
%getClusterCentroids Computes cluster centers based on data & index
% Inputs:
% x - m x n matrix with m observations and n variables
% idx - m x 1 vector indicating the cluster for each observations
% statType - optional input specifying whether the centroid is mean(default)
%           or median
% Output: 
% xCent - y x n matrix with y clusters and n variables
% uIdx - y x 1 vector mapping centroid to corresponding cluster

if(nargin == 2)
    statType = 'mean';
end
[~ , n] = size(x);
uIdx = unique(idx);
xCent = zeros(numel(uIdx),n);
f = zeros(numel(uIdx),1);
for iCls = 1:numel(uIdx)
    ii = idx == uIdx(iCls);
    if(strcmpi(statType,'median'))
        xCent(iCls,:) = median(x(ii,:));
    else
        xCent(iCls,:) = mean(x(ii,:));        
    end
    f(iCls) = sum(ii);
end


end

