function [ mergeVal ] = testClusterMerging( clusterInfo1,clusterInfo2,numColumns,alpha )
%testClusterMerging Summary of this function goes here
%   Detailed explanation goes here


% Get data for cluster 1

numRows = sum(cell2mat(clusterInfo1(:,2)));

d1 = zeros(numRows,numColumns);
cnt = 0;
for i = 1:size(clusterInfo1,1)
    try
    fid = fopen(clusterInfo1{i,1},'r');
%     data = fread(fid,[clusterInfo1{i,2} numColumns],'double');
    d1(cnt+1:cnt+clusterInfo1{i,2},:) = fread(fid,[clusterInfo1{i,2} numColumns],'double');
    fclose(fid);
    catch e
         fclose(fid);
         rethrow(e);
    end
    cnt = cnt+clusterInfo1{i,2};
%     clusterInfo1{i,2}
end
numRows = sum(cell2mat(clusterInfo2(:,2)));
d2 = zeros(numRows,numColumns);
cnt = 0;
for i = 1:size(clusterInfo2,1)
    try
    fid = fopen(clusterInfo2{i,1},'r');
%     data = fread(fid,[clusterInfo1{i,2} numColumns],'double');
    d2(cnt+1:cnt+clusterInfo2{i,2},:) = fread(fid,[clusterInfo2{i,2} numColumns],'double');
    fclose(fid);
    catch e
         fclose(fid);
         rethrow(e);
    end
    cnt = cnt+clusterInfo2{i,2};
%     clusterInfo1{i,2}
end

m1 = median(d1);
m2 = median(d2);

a1 = pdist2(d1,m1,'Euclidean');
a2 = pdist2(d2,m1,'Euclidean');
pp1 = ranksum(a1,a2,'alpha',alpha);

a1 = pdist2(d1,m2,'Euclidean');
a2 = pdist2(d2,m2,'Euclidean');
pp2 = ranksum(a1,a2,'alpha',alpha);

mergeVal = ~and(pp1,pp2);


end

