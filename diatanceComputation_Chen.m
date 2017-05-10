% Distance computation for chen
% Created by Santosh Hariharan
% *********** May 4 2016 *************
clear;clc;
% select root directory
[fileNames,rootDir] = uigetfile('*.txt','Select files for distace computation',...
                        'Multiselect','On');
opFileNamePrefix = 'euclideanDis_';

% Get Header
fileNames = fileNames';
fileNames = fileNames(randperm(length(fileNames)));

try
    fid = fopen(fullfile(rootDir,fileNames{1,:}),'r');
    header = regexp(strtrim(fgetl(fid)),'\t','split');    
    fclose(fid);
%     formatStr = [repmat('%s\t',1,txtEnd) repmat('%f\t',1,numel(hdr) - txtEnd)];
%     fid = fopen(fullfile(rootDir,fileNames),'r');
%     t = textscan(fid,formatStr,'Headerlines',1,'delimiter','\t');
%     fclose(fid);
    
catch expc
    fclose(fid);
    rethrow(expc);
end

[selection, ~] = listdlg('Liststring',header,'PromptString','Select Metadata',...
                        'SelectionMode','Multiple');
textHeader = false(1,numel(header));
textHeader(1,selection) = true;
formatStr = [];
for i = 1:numel(header)
    if(sum(selection == i)>0)
        formatStr = [formatStr '%s '];
    else
        formatStr = [formatStr '%f '];
    end
end
% formatStr = formatStr';

[selection2, ok] = listdlg('Liststring',header(1,textHeader),'PromptString',...
                        'Column for treatments',...
                        'SelectionMode','Single');

% Read Files
columnForGrouping = selection(selection2);
maxrows = 1000000;
textdata = cell(maxrows,1);
data = nan(maxrows,sum(~textHeader));
cnt = 1;
h = waitbar(0,'Reading files....');
for iFiles = 1:numel(fileNames)
    try
        fid = fopen(fullfile(rootDir,fileNames{iFiles,:}),'r');
        fprintf('%s\n',fileNames{iFiles,:});
        t = textscan(fid,formatStr,'Headerlines',1,'delimiter','\t');
        fclose(fid);
    catch expc
        fclose(fid);
        fprintf('Skipped %s\n',fileNames{iFiles,:});
    end
    m = size(t{1,1},1);
    data(cnt:cnt+m-1,:) = cell2mat(t(1,~textHeader));
    textdata(cnt:cnt+m-1,:) = t{1,columnForGrouping};
    cnt = cnt+m;
    waitbar(iFiles./numel(fileNames),h);
end
close(h);
ii = sum(isnan(data),2) == 0;
data = data(ii,:);
textdata = textdata(ii,:);
data = zscore(data);
clear ii cnt m t maxrows textHeader header fileNames
clear selection selection2
%% Compute euclidean distance of means
uTrt = unique(textdata);
[m,n] = size(data);
centroidTrt = nan(numel(uTrt),n);
for i = 1:numel(uTrt)
    ii = strcmpi(uTrt{i,:},textdata);
    centroidTrt(i,:) = mean(data(ii,:));
end

distanceMat = pdist2(centroidTrt,centroidTrt);
m = numel(uTrt);
% Write output file


opFileName = [opFileNamePrefix datestr(now,'yyyymmddHHMMSS') '.txt'];
try
    fid = fopen(fullfile(rootDir,opFileName),'w');
    fprintf(fid,'%s\t','');
    for i = 1:m
        fprintf(fid,'%s\t',uTrt{i,:});
    end
    fprintf(fid,'\n');
    for i = 1:m
        for j = 1:m+1
            if(j==1)
                fprintf(fid,'%s\t',uTrt{i,:});
            else
                fprintf(fid,'%f\t',distanceMat(i,j-1));
            end
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
catch expc
    fclose(fid);
    rethrow(expc)
end
msgbox('Completed analysis');




                    
                    
