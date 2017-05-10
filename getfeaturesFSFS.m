function [ f ] = getfeaturesFSFS( miLambda,k )
%getfeaturesFSFS Selects features based on algorithm by Mitra P. amd Pal S.
% "IEEE Transactions on Pattern Analysis and Machine Intelligence" March
% 2002
% x - mxd matrix of m observations and d variables or features
% Output
% f - indices of features selected

% Check for empty marix
if(isempty(miLambda))
    f = 0;
    return;
end
[~, d] = size(miLambda);
if(d == 1)
    f = 1;
    return;
elseif(k>d || k ==1)
    f=1:d;
    return;
% else
    
end

for i = 1:d
    miLambda(i,i) = inf;
end
% origIX = repmat(1:d,d,1);
[sortedmiLambda,origSort] = sort(miLambda,2);
% sortedmiLambda = sortedmiLambda(:,2:end); % Remove the first column
% origSort = origSort(:,2:end);
f = true(d,1);allfeatures = false(d,1);
% Start looping till all features are selected or rejected
while (d-sum(allfeatures))>0 %Outer iteration
%     Find the first k nearest neighbors for each feature

    if(k<=1)
        break;
    end
    [sortedmiLambda,sortedIx] = sort(sortedmiLambda,2);  
    disp('');
    for ifeat = 1:size(origSort,1)
        sortedIx(ifeat,:) = origSort(ifeat,sortedIx(ifeat,:));
    end
    k = k-1;
    fprintf('Initial K :%i\n',k);
    [rik,ftokeep] = min(sortedmiLambda(:,k),[],1);
    fprintf('Before rik :%f\n',rik);
    err=rik;rik = inf;
    
    f(sortedIx(ftokeep,1:k)) = false;
    allfeatures(sortedIx(ftokeep,1:k)) = true;
    allfeatures(ftokeep,1) = true;
    
    if(k > sum(f)-1) 
        k = sum(f) - 1;
%     else
%         k = k-1;
    end
    sortedmiLambda(ismember(sortedIx,find(~f))) = inf;
    sortedmiLambda(ftokeep,:) = inf;
%     sortedmiLambda(:,ftokeep) = inf;
    sortedmiLambda(sortedIx(ftokeep,1:k),:)  = inf;
%     sortedmiLambda(:,sortedIx(ftokeep,1:k))  = inf;
    fprintf('sum(f) : %d  sum(allfeatures) :%d\n',sum(f),sum(allfeatures));
    while rik>err
        if(k==1 || k==0)
            break;
        end
        [rik,ftokeep] = min(sortedmiLambda(:,k),[],1);  
        fprintf('rik :%f err :%f k :%d \n',rik,err,k);
        if(rik<err)
            fprintf('Entered While');
            f(sortedIx(ftokeep,1:k)) = false;
            sortedmiLambda(ismember(sortedIx,find(~f))) = inf;
            sortedmiLambda(ftokeep,:) = inf;
%             sortedmiLambda(:,ftokeep) = inf; 
            sortedmiLambda(sortedIx(ftokeep,1:k),:)  = inf;
%             sortedmiLambda(:,sortedIx(ftokeep,1:k))  = inf;
            allfeatures(sortedIx(ftokeep,1:k)) = true;
            allfeatures(ftokeep,1) = true;
            fprintf(' sum(f) : %d  sum(allfeatures) :%d\n',sum(f),sum(allfeatures));
        else
            k = k-1;
        end        
    end
    
end
    
f= find(f);    
end












