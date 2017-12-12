function [ node1, node2, edgeWeights,aPicked ] = getGraphByEMST( x,k,cutoff,n )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
if(nargin<1)
    error('@getGraphByEMST: Not enough arguments');
elseif(nargin<2)
    k = size(x,1);n= 500;cutoff = .5;
elseif(nargin<3)
    n= 500;cutoff = .5;
elseif(nargin<4)
    n = 500;
end
if(n)

[aConnections, aPicked] = getWeightMatrixByMST(x,k,n);
aConnections = tril(aConnections,-1);
% jj = aConnections<cutoff;
% aConnections(jj) = 0;
[node1,node2,edgeWeights] = find(aConnections);

end

