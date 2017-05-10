function pval = enrich_score(e_c)
rawmat=e_c;
pval = zeros(size(rawmat));
clustersize = sum(rawmat);
funcsize = sum(rawmat,2);
for i=1:size(rawmat,1)
    for j=1:size(rawmat,2)
        pval(i,j) = 1-hygecdf(rawmat(i,j)-1,sum(funcsize),funcsize(i,1),clustersize(1,j));
    end
end
end