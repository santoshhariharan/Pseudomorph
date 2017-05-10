function [ac,rc] = runClustering(data,grps,feat)

 
[~,n] = size(data);
% C = clsIn(data(:,feat));
idx = kmeans(data(:,feat),numel(unique(grps)),'Replicates',500);
% idx = apclusterK(C.S,numel(unique(grps)));
uIdx =  unique(idx);
for i = 1:numel(uIdx)
    idx(idx==uIdx(i)) = uIdx(i)+103;
end
uIdx =  unique(idx);
for i = 1:numel(uIdx)
    idx(idx==uIdx(i)) = mode(grps(idx==uIdx(i)));
end
c = class_metric(confusionmat(grps,idx));
ac = c.accu;
% nmi(grps,idx)
p = randperm(n,sum(feat));
% C = clsIn(data(:,p));
% idx = kmeans(data(:,feat),numel(unique(grps)),'Replicates',500);
idx = kmeans(data(:,p),numel(unique(grps)),'Replicates',500);
% idx = apclusterK(C.S,numel(unique(grps)));
uIdx =  unique(idx);
for i = 1:numel(uIdx)
    idx(idx==uIdx(i)) = uIdx(i)+103;
end
uIdx =  unique(idx);
for i = 1:numel(uIdx)
    idx(idx==uIdx(i)) = mode(grps(idx==uIdx(i)));
end
c = class_metric(confusionmat(grps,idx));
rc = c.accu;



end