function idx = getFocusFilteredData(data,param)
% Performs filtering based on previously defined RF classifier
% Inputs: Data to be filtered
% data: Data to be filtered
% param: A strcuture that contains meta information and classifier

data = bsxfun(@minus,data, param.focusmean);
data = bsxfun(@rdivide,data,sqrt(param.focusstd));
data = data(:,logical(param.fredfeat));
Y = classRF_predict(data,param.fnet);
idx = false(size(data,1),1);
idx(Y==1) = true;

end