function [ xOut ] = getNonOverlapCoordinates( x, centers, radius )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

k = 10;
range = min(abs(max(x) - min(x)));
addDis = .01*range;
[idx, dis] = knnsearch(x,centers,'K',k+1);
oriX = x;
idx = idx(:,2:end);
dis = dis(:,2:end);

for i = 1:size(centers,1)
    ii = dis(i,:)<=radius(i);
    coord2move = idx(i,ii);
    origDis = dis(i,ii);
    if(~isempty(coord2move))
        for j = 1:numel(coord2move)
            xOld = x(coord2move(j),:);
            v = xOld - centers(i,:);            
            u = v./sqrt(sum(v.^2));
            nD = radius(i)+addDis;
            xNew = centers(i,:) + nD.*u;
            x(coord2move(j),:) = xNew;
        end
    end
    
end
xOut = x;

end

