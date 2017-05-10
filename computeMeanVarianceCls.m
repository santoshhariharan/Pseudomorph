function [md] = computeMeanVarianceCls(pth,md)
% Incremental evaluation of mean and variance of data files.
% The output is saved as meanvar.mat 
% Output Variables
% meaninc - Mean of the data so far
% varinc - variance of the data (NOT STANDARD DEVIATION)
% numpinc - number of rows processed so far
% Computattion of Mean and variace requires a certain folder structure
% Root Folder --> Individual file folders --> Files
% Modified to use textscan as opposed to importdata
% Modified: May 10 2016

% textfeat = md.textfeat;
% datafeat = md.datafeat;
listing = dir(pth);
% meaninc = 0;varinc = 0; numpinc=0;
firstFilenum=false;
colNan = zeros(1,sum(~md.textfeat));
maxColVal = zeros(1,sum(~md.textfeat));
minColVal = zeros(1,sum(~md.textfeat));
% h2 = wb(handles.statusbarpanel);
for kk=3:size(listing)
    if(listing(kk).isdir)
        continue;
    end  
        
    mstart = regexpi(listing(kk).name,'.txt');
    if(isempty(mstart))
        continue;
    end
    
    fprintf('File name %s\n',listing(kk).name);
%     Data = importdata(listing(kk).name);
%     hdr = Data.textdata(1,:);
% Open file 
    fid  = fopen(fullfile(pth,listing(kk).name),'r');
    data = textscan(fid,md.formatString,'headerlines',1,'delimiter','\t');
    fclose(fid);
    data = cell2mat(data(:,~md.textfeat));
    if(isempty(data))
        fprintf('   No data found in file\n');
    end
    infRows = sum(isinf(data),2)>0;
    data(infRows,:)=[];
    nanRows = sum(isnan(data),2)>0;
    colNan = colNan + sum(isnan(data));
    data(nanRows,:)=[]; % Can comment it later if one columns has all Nan's
    
%     Data.textdata(infRows,:)=[];
    
    if(~firstFilenum)
        [md.meaninc, md.varinc, md.numpinc ] = batch_meanvar(data);        
        firstFilenum = true;
%         txthdr = hdr(1,textfeat);
%         datahdr = hdr(1,~textfeat);
    else
        [md.meaninc, md.varinc, md.numpinc ] = batch_meanvar(data,md.meaninc,md.varinc,md.numpinc);
    end
    maxColVal = nanmax(maxColVal,nanmax(data,[],1));
    minColVal = nanmax(minColVal,nanmin(data,[],1));
end
% md.meaninc = meaninc;
% md.varinc = varinc;
% md.numpinc = numpinc;
md.numNanCol = colNan./md.numpinc;
md.maxColVal = maxColVal;
md.minColVal = minColVal;
end
    