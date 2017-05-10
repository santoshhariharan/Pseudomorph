function finalFeatures = unsupervisedFeatureReduction(fileMap,cntrl,param,reductionType)

% perfrom feature reduction on data
% Load 10% of the entire data
% Split randomly into 5 parts
% For each part perform feature reduction
% Pick features that are atleast represented twice
%% Choose Method of feature selection
if(nargin==3)
    reductionType = 'GreedyFS';% Other options are PCA, fsfs
end
%% Feature selection per control and feature selection for entire data
% Pick a maximum of 10000 points per control
numRepeates = 10;
numFeatures = 30;
finalFeatures = false(1,sum(param.datafeat));
tmpFeat = zeros(1,sum(param.datafeat));
maxPoints = 100000;
% data = -1*ones(maxPoints*size(cntrl,1),sum(param.datafeat));
data = zeros(maxPoints,sum(param.datafeat));
% groups = zeros(maxPoints*size(cntrl,1),1);
cnt = 1;
fprintf('Feature selection...............................\n')
for jRepeat = 1:numRepeates
    fprintf('\b\b%2.0f',jRepeat);
    for iControls = 1:size(cntrl,1)
        files = fileMap(strcmpi(fileMap(:,2),cntrl{iControls,:}),1);
        d = readfilesMod(files,param); 
        ind = randperm(size(d,1),round(.6*size(d,1)));
        d = bsxfun(@minus,d, param.meaninc);
        d = bsxfun(@rdivide,d,sqrt(param.varinc));
        data(cnt:cnt+numel(ind)-1,:) = d(ind,param.datafeat);
        cnt = cnt +numel(ind);
    end
%     cnt = cnt -size(D.data,1);
    clear D;
    if(cnt < maxPoints)
        data = data(1:cnt,:);
    end 
    if(strcmpi(reductionType,'pcaselect'))
        redFeatures = unsupervisedPCASelect(data,numFeatures);
    else
        redFeatures = unsupervisedGreedyFS(data,numFeatures);
    end
    finalFeatures = or(finalFeatures,redFeatures');
    tmpFeat = tmpFeat+finalFeatures;
end
fprintf('\n');
if(numRepeates>1)
    finalFeatures = tmpFeat>floor(.2*numRepeates);
end
fprintf('selected Features\n');
fprintf('-----------------\n');
p = param.datahdr(1,param.datafeat);
I = find(finalFeatures);
for i = 1:sum(finalFeatures)
    fprintf('%s\n',p{1,I(i)});
end
clear numRepeates data cnt redFeatures maxPoints files;
clear iControls jRepeat


