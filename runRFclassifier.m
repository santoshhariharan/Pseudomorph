function [ac,rc] = runRFclassifier(data,groups,feat)
    
    numruns = 5;
    [m,n] = size(data);
    ac = zeros(numruns,1);
    rc = zeros(numruns,1);
    uGrp = unique(groups);
    sumfeat = sum(feat);
    % Pick random 1000 points from each group
    for iRuns = 1:numruns
        ind = zeros(1000*numel(uGrp),1);
        cnt = 0;
        for ig = 1:numel(uGrp)
            ii = find(groups == uGrp(ig));
            ind(cnt+1:cnt+1000) = ii(randperm(numel(ii),1000));
            cnt = cnt+1000;
        end
        ii=~ismember(1:m,ind); 
    %     traindata = data(ind,:);
    %     B = TreeBagger(50,data(ind,feat),groups(ind,:),'NvartoSample',sqrt(sum(feat))) ;
        B = classRF_train(data(ind,feat),groups(ind,:),50,sqrt(sumfeat));           
        ghat = classRF_predict(data(ii,feat),B);
        c = class_metric(confusionmat(groups(ii,1),ghat));
        ac(iRuns) = c.accu;
        
%     Random features
        ix = randperm(n,sumfeat);
        B = classRF_train(data(ind,ix),groups(ind,:),50,sqrt(sumfeat));      
        ghat = classRF_predict(data(ii,ix),B);
        c = class_metric(confusionmat(groups(ii,1),ghat));
        rc(iRuns) = c.accu;

    end
    
end