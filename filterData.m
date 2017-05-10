function [ indexToRemove ] = filterData( data, param )
%filterData Performs various data specific filtering & retuns cleaned data
%   Detailed explanation goes here


if(size(data,2) ~= numel(param.datahdr))
    indexToRemove = false(size(data,1),1);
    disp('@filterData: No filtering performed');
    return;
end
indexToRemove = false(size(data,1),1);


% Focus filtering
if(isfield(param,'fnet'))
    if(~isempty(param.fnet))        
        indexToRemove = or(indexToRemove,~getFocusFilteredData(data,param));
    end
end

% Intensity Filtering
if(isfield(param,'intensityLowerThreshold'))
    indexToRemove = or(indexToRemove,or(data(:,param.intensityChannel) <= param.intensityLowerThreshold,...
        data(:,param.intensityChannel) >= param.intensityUpperThreshold));
end

% Morphometric Filtering
if(isfield(param,'morphometricFilterLowerThreshold'))
    d = data(:,param.morphChannel);
    d = d(:,2)./d(:,1);
    indexToRemove = or(indexToRemove,or(d <= param.morphometricFilterLowerThreshold,...
        d >= param.morphometricFilterUpperThreshold));
end

end

