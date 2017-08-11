function [ A, W ] = createGraph(index,dis)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


m = size(index,1);
A = zeros(m);
W = zeros(m);

for i = 1:m
    A(i,index(i,:)) = 1;
    W(i,index(i,:)) = dis(i,:);
end


% Create symmetric Adj matrix
A = A+A';
A = double(A>0);

W = W+W';
W = W/2;
end

