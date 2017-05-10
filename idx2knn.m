function  G  = idx2knn(IDX,D)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

sigma = std(D(:));
D = exp((-D.^2)/(2*sigma*sigma));
fprintf('@idx2knn\n')
[n,k] = size(IDX);
I = nan(1, k*n );
J = I;
S = I;
row = 1;
pctdone = 10;
ticker = round(n/10);
t = tic;
for ii = 1:n
    pt_neighbs = IDX(ii,:);
    idx = row:row+k-1;
    I(idx) = repmat(ii,1,k);
    row = row+k;
    J(idx) = pt_neighbs;
    S(idx) = D(ii,:);
    if ~mod(ii,ticker) && ii>1
        fprintf(1,'%i percent complete: %.2f s\n',...
            pctdone, toc(t) );
        pctdone = pctdone + 10;
    end
end
G = sparse(I,J,S,n,n);
clear I J S IDX
u = triu(G,1);
v = tril(G,-1);
G = v + u';
end

