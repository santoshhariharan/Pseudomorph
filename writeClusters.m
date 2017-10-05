function  writeClusters( ind,text,filename,wFlag,mode)
% Outputs the percentage of Clusters on the screen.
% Plot enables the cluster percentages to be plot in a Bar graph
% Inputs:
% ind = mx1-m: number of data points. Values indicate the index of the
%       centroids
% grps= mx1: Group numbers of each point
% data = mxd array-d is hte dimensions of the array
% text - Text data for output
% PerCluster - Output type. If ture indicates the output per cluster. Else indicates
%       output per Sample
% pcomp - Turn on principal component plot
% Gather Inputs
if(nargin == 3)
    mode = 1;
    wFlag = 'overwrite';
elseif(nargin ==4)
    mode = 1;
elseif(strcmpi(mode,'separate'))
    mode = 0;
else
    mode =1;
    wFlag = 'overwrite';
end


numClusters = unique(ind);

oldF = pwd;
if(exist(['output_' date],'dir') ~= 7)
    mkdir([pwd '\output_' date]);
%     disp('DDDD')
end
cd([pwd '\output_' date]);
    
% end
% mkdir(['output_' date]);
% cd('./Output')
% if(exist(['output_' date],'dir') ~= 7)
%     mkdir(['output_' date]);
% end
% cd(['./output_' date]);
c1 = {'Well-No.' 'Plate-ID' 'Image-No.' 'Control' 'Field-of-View' 'X-Coord' 'Y-Coord' 'Classes'};
% clas = zeros(size(ind,1),1);
clas = zeros(numel(ind),1);
text = text(:,[1 2 3 9 12 13 14]);temp={};
for i=1:numel(numClusters)
    f = find(ind == numClusters(i));    
    clas(f) = i;      
%     f(I)
    if(~mode)
        temp = [text(f,:) strtrim(cellstr(num2str(clas(f))))];
        temp = [c1;temp];
       
        writestr([filename '_Cluster' num2str(i) '_output_' date '.txt'],temp,'Append');
    else
%         temp = [temp;text(f,:)];
    end
end 

if(mode)
    temp = [temp;text];
    text = [temp strtrim(cellstr(num2str(clas)))];
   
    if(exist([filename '_output_' date '.txt'],'file'))
        if(strcmpi(wFlag,'Append'))
            c2 = [text];
        else
            c2 = [c1;text];
        end
    else
        c2 = [c1;text];
    end
%     c1 = [c1;text];
    writestr([filename '_output_' date '.txt'],c2,wFlag);
end

cd(oldF);
end

