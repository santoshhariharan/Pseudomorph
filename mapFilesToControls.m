function fileMap = mapFilesToControls(pth,cntrl)
% Map individual files to controls
% Requires user to be in the corrrect directory
% currDir = pwd;
l = dir(pth);
% cntrl = {'ergic';'golgin';'tabik'};
fileMap = cell(size(l,1),1);
tmp = false(size(l,1),1);
for iFiles = 3:size(l,1)
    if(l(iFiles).isdir)
        continue;
    end    
    mstart = regexpi(l(iFiles).name,'.txt');
    if(isempty(mstart))
        continue;
    end
    token = strsplit(l(iFiles).name,'_');
    for jCnt = 1:numel(cntrl)
 
        if(sum(strcmpi(token,cntrl{jCnt,:}))>0)
            fileMap{iFiles,1} = l(iFiles).name;
            fileMap{iFiles,2} = cntrl{jCnt,:};
            tmp(iFiles) = true;
        end
    end
end
fileMap = fileMap(tmp,:);
clear tmp iFiles jCnt l ii currDir mstart token
end