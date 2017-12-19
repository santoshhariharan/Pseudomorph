%% Mean and Deviation compute
% 
clc;clear all;
 fprintf('Starting mean computation\n');
% Define Inputs:
% pth='F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data';
pth = 'F:\Projects\Proteinlocalization\PseudoMorph\Bin2Data\AllControlFiles';


columnForControls = 9;
fprintf('Path: %s\n',pth);
fprintf('Control column: %i\n',columnForControls);
load(fullfile(pth,'parameters.mat'));% Load parameter file
param.rootpath = pth;
intensityFeature = 'Ch2_INT_Cell_intensity';
intFeat = strcmpi(intensityFeature,param.datahdr);
nucAreaFeat = strcmpi('Ch1_MOR_Nucleus_area',param.datahdr);
cellAreaFeat = strcmpi('Ch1_MOR_Cytoplasm_area',param.datahdr);
filePrefix = '.txt';
fNames = dir(pth);
% Module 3: Read & Load data after filtering
fprintf('Module 3.......\n');
mxRw = 1000000; 
allD = zeros(mxRw,sum(param.datafeat));
allInten = zeros(mxRw,1);
allMorRatio = zeros(mxRw,1);
allTxt = cell(mxRw,1);


cnt = 1;
fprintf('Completed Reading................');
for iFiles = 3:size(fNames,1)
    fprintf('\b\b\b\b\b\b\b\b\b%8.3f%%',iFiles*100./size(fNames,1));    
    if(fNames(iFiles).isdir)
        continue;
    end
    tok = regexpi(fNames(iFiles).name,filePrefix,'match');
    if(isempty(tok))
        continue;
    else
        D = readfiles(cellstr(fNames(iFiles).name),param);
    end
    %     Remove cells out of focus
    focus = getFocusFilteredData(D.data,param);
    D.data = D.data(focus,:);
    D.textdata = D.textdata(focus,:);
    ii = (D.data(:,strcmpi('Ch1_INT_Nucleus_intensity',param.datahdr))./...
        D.data(:,strcmpi('Ch1_INT_Cytoplasm_intensity',param.datahdr)))>3.5;
    jj = (D.data(:,strcmpi('Ch1_INT_Nucleus_intensity_stddev',param.datahdr))./...
        D.data(:,strcmpi('Ch1_INT_Cytoplasm_intensity_stddev',param.datahdr)))>3.5;
    ii = and(ii,jj);  
    
    D.data = D.data(ii,:);
    D.textdata = D.textdata(ii,:);
    allInten(cnt:cnt+size(D.data,1)-1,:) = D.data(:,intFeat); 
%     allMorIntensity(cnt:cnt+size(D.data,1)-1,:) = D.data(:,roiIntFeat);
    allMorRatio(cnt:cnt+size(D.data,1)-1,:) = D.data(:,nucAreaFeat)./...
                                            (D.data(:,cellAreaFeat)+D.data(:,nucAreaFeat));    
    allD(cnt:cnt+size(D.data,1)-1,:) = D.data(:,param.datafeat);
    allTxt(cnt:cnt+size(D.data,1)-1,:)= D.textdata(:,columnForControls);
    cnt = cnt+size(D.data,1);
end
if(cnt<mxRw)
    allD = allD(1:cnt-1,:);
    allInten = allInten(1:cnt-1,:);
    allTxt = allTxt(1:cnt-1,:);
    allMorRatio = allMorRatio(1:cnt-1,:);
end
fprintf('\n');
%% Remove Artifacts and noise

% Remove NaN entries
ii = sum(isnan(allD),2) ==0;
allD = allD(ii,:);
allInten = allInten(ii,:);
allTxt = allTxt(ii,:);
allMorRatio = allMorRatio(ii,:);
% allTxtValues = allTxtValues(ii,:);
% allMorIntensity = allMorIntensity(ii,:);
% allTxtOrg= allTxtOrg(ii,:);
fprintf('\nRemoved NAN %i\n',sum(~ii));

%% Remove incorrectly segmented & low intensity Objects
ii = allMorRatio <= .5 & allMorRatio >= .2; % Value of this needs optimization
allD = allD(ii,:);
allTxt = allTxt(ii,:);
allInten = allInten(ii,:);
% allTxtOrg = allTxtOrg(ii,:);
% allInten = allInten(ii,:);
% allMorRatio = allMorRatio(ii,:);
% allTxtValues = allTxtValues(ii,:);
% allMorRatio = allMorRatio(ii,:);
% allMorIntensity = allMorIntensity(ii,:);
fprintf('#Cells Removed 4 Morphology %i of %i\n',sum(~ii),numel(ii));


%% Remove cells with low intensity
ii = allInten > 100;
allD = allD(ii,:);
allTxt = allTxt(ii,:);
% allTxtOrg = allTxtOrg(ii,:);
allInten = allInten(ii,:);
% allTxtValues = allTxtValues(ii,:);
fprintf('#Cells Removed 4 Intensity %i of %\n',sum(~ii));

% Remove Lower 5% and upper 5% data for each control
uControls = unique(allTxt);
jj = false(size(allTxt,1),1);
for i = 1:numel(uControls)
    ii = find(strcmpi(uControls{i,:},allTxt));
    kk = allInten(ii,1)>quantile(allInten(ii,1),.05) &...
                    allInten(ii,1)<quantile(allInten(ii,:),.95);
    jj(ii(kk)) = true;
end
allD = allD(jj,:);
allTxt = allTxt(jj,:);
allInten = allInten(jj,:);
% param.meaninc = mean(allD);
% param.varinc = var(allD);
fprintf('#Cells removed by lower-upper quartile %i\n',sum(~jj));


% % Remove intensity correlated features
% rho = corr(allD,allInten);
% newHeader = param.datahdr(1,param.datafeat);
% ii = find(param.datafeat);
% corrFeat = rho>-.5 & rho < .5;% Retain columns between -.5 & 0,5
% param.datafeat(1,ii(~corrFeat)) = false;
% allD = allD(:,corrFeat);
% fprintf('#Features removed due to correlation %i\n',sum(~corrFeat));
% fprintf('%s\n',newHeader{1,~corrFeat});
% newHeader = newHeader(1,corrFeat);


% % Remove features having 75% same data
% [~,F] = mode(allD,1);
% F= F./size(allD,1);
% jj = F<.75;
% ii = find(param.datafeat);
% param.datafeat(1,ii(~jj)) = false;

% newHeader = newHeader(1,jj);
% allD = allD(:,jj);
% Feature Selection/Reduction
% newHeader = param.datahdr(1,param.datafeat);

% if(featureReduction)    
% %     redFeatures = unsupervisedGreedyFS(allD,numFeatures);
%     [redFeatures,scores] = unsupervisedPCASelect( allD,numFeatures );
% else
%     redFeatures = true(1,size(allD,2));
% end

% [cf]= princomp(allD(:,redFeatures));
% allD = allD(:,redFeatures);
% Print number of cells per control
for i = 1:numel(uControls)
    ii = (strcmpi(uControls{i,:},allTxt));
    fprintf('%s\t: %d\n',uControls{i,:},sum(ii));    
end

% Normalization
param.meanD = mean(allD);
param.stdD = std(allD);
allD = zscore(allD);  

% Save Mean and STD


clear D focus cnt tok iFiles mxRw 
clear allMorRatio 
clear ii jj kk rho mxRw 
clear intensityFeature roiIntFeat
clear fNames filePrefix randPrc
