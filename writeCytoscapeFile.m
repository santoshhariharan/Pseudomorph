function [ output_args ] = writeCytoscapeNetworkFile(netData,nodeData,filename )
%writeCytoscapeFile Write cytoscape files


m = size(nodeData,1);

try
    fid = fopen(filename,'w');
    for i = 1:size(netData,1)
        fprintf('%s\t%s\t',['N' num2str(netData(i,1))],['N' num2str(netData(i,2))]);
        fprintf('%f\n',netData(i,3));
    end    
    fclose(fid);
catch expc
    fclose(fid);
    rethrow(expc)
end


end

