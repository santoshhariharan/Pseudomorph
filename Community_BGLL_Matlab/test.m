clc;clear all;
m = randn(100000,100);
tic;
% idx = knnsearch(m,m,'k',21);
% idx  = idx(:,2:end);
% M = knn2jaccard(idx);
% cluster_jl_cpp(full(M));
phenograph(m,20);
toc;