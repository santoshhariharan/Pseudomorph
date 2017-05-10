function [ h,p ] = mergeClusterByPermutation( c1,c2,alpha,mergeType )
%mergeClusterByPermutation Tests whether two clusters could be merged by
%permutation
% Inputs:
%   c1  - m1 x n data from cluster 1 with m1 observations and n variables
%   c2 - m2 x n data from cluster 2 with m2 observations and n variables
% Output:
%   h - 0 if the clusters cannot be merged and 1 if the clusters can be
%   merged
%   p - P value for merge test
% if(nargin == 2)
%     alpha = .05;
%     mergeType = 'mean';
% end
numPermutations = 500;
h = 1;
% p = 1;
[m1, ~] = size(c1);
[m2] = size(c2,1);
dataMat = [c1;c2];
if(strcmpi(mergeType,'mean'))
    n = m1+m2;
    mOrig = (m1./n)*mean(pdist2(c1,mean(c1))) + (m2./n)*mean(pdist2(c2,mean(c2)));
    mPerm = zeros(numPermutations,1);
    for i = 1:numPermutations
        grps = assignRandomGroups(m1+m2,m1);        
        mPerm(i) = (m1./n)*mean(pdist2(dataMat(grps==1,:),mean(dataMat(grps==1,:)))) +...
                    (m2./n)*mean(pdist2(dataMat(grps==0,:),mean(dataMat(grps==0,:))));
    end
    p = sum(mPerm<mOrig)./numPermutations;
else
    % Compute the original within cluster variance & between cluster
    xminusmean = bsxfun(@minus,c1,mean(c1));
    yminusmean = bsxfun(@minus,c2,mean(c2));
    
    wOrig = (1/m1).*xminusmean'*xminusmean + (1/m2).*yminusmean'*yminusmean;
    wOrig = det(wOrig);
    
    tic;
    % Run number of permutations and create a vextor of w
    wPerm = zeros(numPermutations,1);
    % grps = uint8(zeros(m1+m2,1));
    for i = 1:numPermutations
        grps = assignRandomGroups(m1+m2,m1);
        xminusmean = bsxfun(@minus,dataMat(grps==1,:),mean(dataMat(grps==1,:)));
        yminusmean = bsxfun(@minus,dataMat(grps==0,:),mean(dataMat(grps==0,:)));
        wPerm(i) = det((1/m1).*xminusmean'*xminusmean + (1/m2).*yminusmean'*yminusmean);
    end
    p = sum(wPerm<wOrig)./numPermutations;
end
if(p<=alpha)
    h = 0;
end
fprintf('@mergeClusterByPermutation: Time taken %.3fs\n',toc);
end

function g = assignRandomGroups(n,n1)
g = zeros(n,1);
g(randperm(n,n1)) = 1;
end

