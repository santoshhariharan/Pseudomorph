function [ c ] = computeDistribution( index,labels )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

uI = unique(index);
uL = unique(labels);

c = zeros(numel(uI),numel(uL));

for i = 1:numel(uI)
    ii = index==uI(i);
    for j = 1:numel(uL)
        jj = strcmpi(labels,uL{j,:});
        c(i,j) = sum(ii.*jj);
    end
end
end

