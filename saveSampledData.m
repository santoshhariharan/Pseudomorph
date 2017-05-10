% Get random data from controls ~1000 points from each control
% Make PCA plot and view it
% Save data
pth = 'F:\Projects\Proteinlocalization\PseudoMorph\TestDataForPseudoMorph';
cntrl = {'golgin';'ergic';'tabik'};
fileMap = mapFilesToControls(pth,cntrl);
load(fullfile(pth,'parameters.mat'));

% Run through the controls and write out data
% Reduce data and save them as mat file
% Save as minimum data as possible
cnt = 1;
pointsToSelect = 1500;
maxrows = 4000;
data = zeros(maxrows,sum(param.datafeat));
textdata = cell(maxrows,1);
for iFiles = 1:numel(cntrl)
    ii = strcmpi(fileMap(:,2),cntrl{iFiles,:});
    [d,t] = readfilesMod(fileMap(ii,1),param);
    n = randperm(size(d,1),pointsToSelect);
    d = bsxfun(@minus,d, param.meaninc);
    d = bsxfun(@rdivide,d,sqrt(param.varinc));
    data(cnt:cnt+numel(n)-1,:) = d(n,param.datafeat);
    textdata(cnt:cnt+numel(n)-1,:) = t(n,9);
    cnt = cnt+numel(n);
end
clear d t;
if(cnt<maxrows)
    data = data(1:cnt,:);
    textdata = textdata(1:cnt,:);
end
save(fullfile(pth,'sampleDataCellFeat.mat'),'data','textdata');
return;
% [~,score] = princomp(data);
score = compute_mapping(data, 't-SNE', 3);
score = score(:,1:3);


%%
h=figure; hold on;
mp = jet(numel(cntrl));
for i = 1:numel(cntrl)
    ii=strcmpi(textdata,cntrl{i,:});
    plot3(score(ii,1),score(ii,2),score(ii,3),'o',...
        'MarkerEdgeColor','None','MarkerFaceColor',mp(i,:));
end
hold off;
xlabel('t-SNE1');ylabel('t-SNE2');zlabel('t-SNE3');
grid on;title('Sampled Data');
legend(cntrl);
savefig(h,fullfile(pth,'sampleDataFiguretSNE.fig'));
close(h);
clear;

