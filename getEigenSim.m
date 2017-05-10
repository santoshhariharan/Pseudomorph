function [lambda] = getEigenSim(x)


varX = var(x);
corrX = corr(x);
% assignin('base','corrX',corrX)
% assignin('base','varX',varX)

lambda = zeros(size(corrX));
n = size(corrX,1);

% Compute lambda
for i = 2:n
    for j=1:i-1
        tmp = varX(i) + varX(j);        
        lambda(i,j) = (tmp) - (sqrt(tmp.^2 - (4.*varX(i).*varX(j).*(1-corrX(i,j).^2))));
    end
end

lambda = lambda+lambda';
lambda = .5.*lambda;

end