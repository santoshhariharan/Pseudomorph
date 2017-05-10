function [ sim ] = getFSFSSim( x )
%getFSFSSim Computes mutual information index for data x
if(isempty(x))
    sim = 0;
    return;
end
[~, d] = size(x);
if(d == 1)
    sim = 0;
    return;
end

% Remove nan rows from x
x = x(sum(isnan(x),2) ==0,:);

% remove completely correlated features
varx = nanvar(x);

varxy = varx'*varx;
rhox = corr(x);

% miLambda = zeros(d);
varxplussq = bsxfun(@plus,dot(varx,varx,1)',dot(varx,varx,1))+2*(varxy);
sim = sqrt(varxplussq) - sqrt(varxplussq - (4.*varxy.*(1-(rhox.^2))));


end

