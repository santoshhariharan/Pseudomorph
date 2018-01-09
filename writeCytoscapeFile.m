function  writeCytoscapeFile(fileprefix,directory,adjacencyMatrix,...
                                                edgeWeights,edgeAttribute,...
                                                clusterIndex,treatmentIndex )
%writeCytoscapeFile Write cytoscape files



networkFilename = [fileprefix '_Networkfile.txt'];
networkFilename = fullfile(directory,networkFilename);

nodeFilename = [fileprefix '_Nodefile.txt'];
nodeFilename = fullfile(directory,nodeFilename);

uIndx = unique(clusterIndex);
uTreatmentIndex = unique(treatmentIndex);
distMat = zeros(numel(uIndx),numel(uTreatmentIndex));

for i = 1:numel(uIndx)
    ii = clusterIndex == uIndx(i);
    for j = 1:numel(uTreatmentIndex)
        jj = treatmentIndex == uTreatmentIndex(j);
        distMat(i,j) = sum(ii.*jj);
    end
end
distMat = bsxfun(@rdivide,distMat,sum(distMat,2));
adjacencyMatrix = tril(adjacencyMatrix);
[r, c] = find(adjacencyMatrix);

% Write Network file
m = numel(r);
try
    fid = fopen(networkFilename,'w');
%     write Header
    fprintf(fid,'%s\t','Source');
    fprintf(fid,'%s\t','Target');
    fprintf(fid,'%s\t','Value');
    fprintf(fid,'%s\t','Weight');
    fprintf(fid,'\n');
    for i = 1:m
        fprintf(fid,'%s\t',num2str(r(i)));
        fprintf(fid,'%s\t',num2str(c(i)));
        fprintf(fid,'%f\t',(edgeWeights(r(i),c(i))));
        fprintf(fid,'%f\t',(edgeAttribute(r(i),c(i))));
        fprintf(fid,'\n');
    end    
    fclose(fid);
catch expc
    fclose(fid);
    rethrow(expc)
end


% Compute Nodes
rNodes = unique(r);
nodeAtributes = zeros(numel(rNodes),max(treatmentIndex));
for i = 1:numel(rNodes)
    ii = find(uIndx == rNodes(i));    
    if(~isempty(ii))
        nodeAtributes(i,:) = distMat(ii,:);
    else
        nodeAtributes(i,treatmentIndex(rNodes(i))) = 1;
    end
end

try
    fid = fopen(nodeFilename,'w');
%     write Header
    fprintf(fid,'%s\t','Name');
    for i = 1:size(nodeAtributes,2)
        fprintf(fid,'%s\t',['TRT' num2str(i)]);
    end    
    fprintf(fid,'\n');    
    for i = 1:numel(rNodes)        
        fprintf(fid,'%s\t',num2str(rNodes(i)));
        for j = 1:size(nodeAtributes,2)
            fprintf(fid,'%f\t',nodeAtributes(i,j));
        end
        fprintf(fid,'\n');
    end    
    fclose(fid);
catch expc
    fclose(fid);
    rethrow(expc)
end

end

