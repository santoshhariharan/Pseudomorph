function [ y ] = unsupervisedGreedyFS( x,k )
%unsupervisedGreedyFS Unsupervised feature selection based on Farahat et
%al. "An Efficient Greedy Method for Unsupervised Feature Selection" 2011
% IEEE,International Conference on Data mining
% Inputs:
% x: m x n matrix with m observations and n variables
% k: Number of output features to be selected
% Output:
% y: Index of features selected

if(nargin < 2)
    y = [];
    disp('@unsupervisedGreedyFS: Not enough arguments');
    return;
end

if(isempty(x))
    y = [];
    disp('@unsupervisedGreedyFS: Empty data matrix');
    return;
end
[~, n] = size(x);
if(k>=n)
    y = [1:n]';
    return;
end

fi = sum((x'*x).^2)';
gi = diag(x'*x);
% Initialize variables
y = false(n,1);omega = zeros(n,k);

aTa = x'*x;
for iFeatures = 1:k
%     fprintf('Iteration: %i\n',iFeatures);
    tmp = fi./gi;
%     tmp(isinf(tmp)) = -inf;
    tmp(y) = -inf; % Forcing features selected to be -inf
    [~,l] = max(tmp);
    y(l) = true;
%     fprintf('#Iteration %i #Feature: %i fi/gi: %f\n',iFeatures,l,tmp(l));
    %     bsxfun(@times,omega(:,1:iFeatures-1),omega(l,1:iFeatures-1));
    if(iFeatures == 1)
        delta_t = x'*x(:,l);
%         gamma_t = B'*x(:,l);
    else
        t1 = bsxfun(@times,omega(:,1:iFeatures-1),omega(l,1:iFeatures-1));
        delta_t =  x'*x(:,l) - sum(t1,2);
%         t1 = bsxfun(@times,psi(:,1:iFeatures-1),omega(l,1:iFeatures-1));
%         gamma_t = B'*x(:,l) - sum(t1,2);
        
    end
    omega(:,iFeatures) = delta_t./sqrt(delta_t(l));
%     psi(:,iFeatures) = gamma_t./sqrt(delta_t(l));
    %  Update fi and gi - Here ifeatures is (t-1) already
    t1 = omega(:,iFeatures).^2;
    gi = gi - t1;
    omegasq = sum(t1).*omega(:,iFeatures).^2;
    t1 = aTa;
    if(iFeatures>1)
        tmp1 = 0;
        for i = 1:iFeatures-1
            tmp1 = tmp1+omega(:,i)*omega(:,i)';
        end
        t1 =  t1 - tmp1;
    end
    t1 = t1*omega(:,iFeatures);
    t1 = 2.*omega(:,iFeatures).*t1;
    fi = fi - t1 + omegasq;
end
% disp('@@@@');
end
