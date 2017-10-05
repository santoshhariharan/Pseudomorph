function [D,err] = readfiles(fnames,param)
% Reads data from tab limited text files
% D.data=[];D.textdata = [];err=0;
if(isempty(fnames))
    D=[];
    err=1;
    error('@readfiles: No files to read');
end
% data = Nan(param.numpinc,sum(~param.textfeat));
% textdata = cell(param.numpinc,sum(param.textfeat));
textend = uint8(sum(param.textfeat));
% fnames = fnames';
D.data = [];D.textdata = {};
% maxData = 0;
for iFiles = 1:size(fnames,1)
%     fnames{i,:}

% Open file
       fid = fopen(fullfile(param.rootpath,fnames{iFiles,:}));
       filedata = textscan(fid,param.formatString,'headerlines',1,'delimiter','\t');
       fclose(fid);
       filetext = {};
       for jText = 1:textend
           filetext = [filetext filedata{:,jText}];
%            textdata(maxData+1:(maxData + size(filedata{:,1},1)),jText) = filedata{:,jText};
       end
       filedata = cell2mat(filedata(:,~param.textfeat));
%        maxData  = maxData + size(filedata{:,1},1);
       D.data = [D.data;filedata] ;
       D.textdata = [D.textdata;filetext];
end
rowsToRemove = sum(isinf(D.data),2)>0 | sum(isnan(D.data),2)>0;
D.data = D.data(~rowsToRemove,:);
D.textdata = D.textdata(~rowsToRemove,:);
err = 0;
end
