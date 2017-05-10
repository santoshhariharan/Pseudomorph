function [ nIndex ] = reindexClusters( ind )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

uIndx = unique(ind);
nIndex = ind;
for i = 1:numel(uIndx);
    nIndex(ind==uIndx(i),1) = i;
end
end

