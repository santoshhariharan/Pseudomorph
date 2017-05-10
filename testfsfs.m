% pp = [1:3];
% vv = bsxfun(@plus,dot(pp,pp,1)',dot(pp,pp,1))+2*(pp'*pp);
clc;close all
% x = randn(300,1);
% x = [rand(300,1) x randn(300,3) 2*x];
load ionosphere;
x = X;
x = x(:,[1 3:end]);
y = unsupervisedGreedyFS(x,20);
%% %
fprintf(1,'Computing Feature Similarities..\n');
dm = zeros(size(x,2));
no_feature = size(x,2);
for i=1:no_feature,
   fprintf(1,'Similarity Computed for Feature No. %d\n',i);
   for j=1:no_feature,
%       x1=data(:,i);x2=data(:,j);
      if i < j
         dm(i,j)=f2f(x(:,i),x(:,j),3);
      else
         dm(i,j)=0.0;
      end
   end
end

% load iris_dataset
% x = irisInputs';
% x = getappdata(0,'alldata');
%%
clc;
[~,n] = size(x);
k = [1:5:n-1];
f= zeros(numel(k),1);
diffFeat = zeros(numel(k),1);

sim = getFSFSSim( x );
% sim = sim(sum(isnan(sim)))
for ik = 1:numel(k)
    tic;
%     t=getfeaturesFSFS(sim,k(ik));
%     t = fsfs(x,no_feature,k(ik));
    t = fsfs_mod(dm,no_feature,k(ik));
    fprintf('Number of features %i Time: %.2fs K: %i\n',numel(t),toc,k(ik));
%     toc;
    f(ik) = numel(t);
    if(ik>1)
        diffFeat(ik) = f(ik) - f(ik-1);
    end
%     re;
end
plot(k(1:end),f(1:end));
% figure;plot(diffFeat);
disp('done');
