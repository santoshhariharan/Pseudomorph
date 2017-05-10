function [ newIndex ] = resetIndexValues( oldIndex )
%resetIndexValues Takes old categories and maps them into 1 to n categories
%Input:
% oldIndex: m x 1 vector with m observations
% Output:
% newIndex: Same as old index but with new categories

newIndex = zeros(size(oldIndex,1),1);
oldIndex = -oldIndex;
uOld = unique(oldIndex);
for i = 1:numel(uOld)
    ii = oldIndex == uOld(i);
    newIndex(ii) = i;    
end

end

