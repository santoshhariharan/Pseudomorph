function writeImageInformation(textInfo,classNum,opFileName)

if(size(textInfo,1) ~= size(classNum,1))
    error('Uneven number of rows');
end

[m,n] = size(textInfo);
try
    fid = fopen(opFileName);
    for i = 1:m        
        for j = 1:n
            fprintf('%s\t',textInfo{i,j});
        end
        fprintf('%f\n',classNum(i));
    end
    fclose(fid);
catch expc
    fclose(fid);
    rethrow(expc);
end


end
