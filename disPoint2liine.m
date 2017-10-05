function [ dis ] = disPoint2liine( x,cord1,cord2 )
%disPoint2liine Computes distance of multiple points to many lines
%   Computes distance of points defined in x to all the lines in L
% x - m x d with m rows and d dimensions
% cord1 - n x d x 2 with n rows and d dimensions and the third dimension 
%       defines the corodinates
% cord2 - n x d with n rows and d dimensions and the third dimension 
%       defines the corodinates
% 
% Ouput:
% 
% dis - m x n - Each value shows the euclidean distance between the
% points and the line



if(size(x,2) ~= size(cord1,2))
    error('@disPoint2liine: Dimension mismatch');
end
[m,~] = size(x);
n = size(cord1,1);

% Distance between points from the line/Length of line
% lineD = L(:,:,1) - L(:,:,2);
% linDSquared = dot(lineD,lineD,2);
dis = inf(m,n);
for i = 1:n
    V1 = repmat(cord1(i,:),m,1);
    V2 = repmat(cord2(i,:),m,1);
    a = V2 - V1;
    b = x - V1;
    aSquared = sqrt(dot(a,a,2));
    bSquared = sqrt(dot(b,b,2));
    abSq = aSquared.*bSquared;
    t  = dot(a,b,2)./(abSq);
    ii = t>0 & t<=1; %Angle between 0 and 90
    theta = acos(t);
    bb = sin(theta).*bSquared;
    dis(ii,i) = bb(ii);
end

end

