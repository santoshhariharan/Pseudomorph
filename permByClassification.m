function [ p,accuOrig,accuRand ] = permByClassification( x1,x2,n )
%permByClassification Perform permutation test based on classification data
%   Detailed explanation goes here

if(nargin == 2)
    n = 100;
end
kk = 5;
m1 = size(x1,1);
m2 = size(x2,1);

minSamples = floor(min(m1,m2));
accuOrig = nan(n,1);
accuRand = nan(n,1);
d = [x1;x2];
lb = [ones(m1,1);2*ones(m2,1)];
% Compute original accuracy
for iRpt = 1:n
    ii = getSamples(lb,minSamples);
    if((sum((lb(ii.training) == 1)) < kk) ||...
            (sum((lb(ii.training) == 2)) < kk))
        continue;
    end
    cl = knnclassify(d(ii.test,:), d(ii.training,:), lb(ii.training),kk);
    cl = confusionmat(lb(ii.test),cl);
    accuOrig(iRpt) = sum(diag(cl))./sum(cl(:));
    
    nLbl = lb(randperm(numel(lb)));    
    ii = getSamples(nLbl,minSamples);
    cl = knnclassify(d(ii.test,:), d(ii.training,:), nLbl(ii.training),kk);
    cl = confusionmat(nLbl(ii.test),cl);
    accuRand(iRpt) = sum(diag(cl))./sum(cl(:));
end

ii = isnan(accuOrig);
accuOrig = accuOrig(~ii);
accuRand = accuRand(~ii);
n = n - sum(ii);
% % Compute permuted accuracy
% for iRpt = 1:n
%     nLbl = lb(randperm(numel(lb)));    
%     ii = getSamples(nLbl,minSamples);
%     cl = knnclassify(d(~ii,:), d(ii,:), nLbl(ii));
%     cl = confusionmat(nLbl(~ii),cl);
%     accuRand(iRpt) = sum(diag(cl))./sum(cl(:));
% end
if(n>1)
%     p = sum((accuRand>=mean(accuOrig)))./n;
    
    [~,p] = ranksum(accuRand,accuOrig);
    accuRand = mean(accuRand);
    accuOrig = mean(accuOrig);
else
    p = 1;
    accuRand = 0;
    accuOrig = 0;
end
end

function sSize = getSamples(grp,minPerGrp)
sSize.training = false(numel(grp),1);
sSize.test = false(numel(grp),1);

uGrp = unique(grp);
trMinPerGrp = floor(.7.*minPerGrp);
tstMinPerGrp = floor(.3.*minPerGrp);
for i = 1:numel(uGrp)
    ii = find(grp ==uGrp(i));
    jj = randperm(numel(ii));
    sSize.training(ii(jj(1:trMinPerGrp))) = true;
    sSize.test(ii(jj(trMinPerGrp+1:trMinPerGrp+tstMinPerGrp))) = true;
%     ii = ii()
%     sSize(ii(jj)) = true;
end


end



