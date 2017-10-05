function [ projectedCoord,relDis ] = projectPointonLine( x,cord1,cord2 )
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

if(size(x,1) ~= size(cord1,1))
    error('@disPoint2liine: Dimension mismatch');
end
[n,~] = size(x);
% n = size(cord1,1);

% Distance between points from the line/Length of line
% lineD = L(:,:,1) - L(:,:,2);
% linDSquared = dot(lineD,lineD,2);
relDis = nan(n,1);
projectedCoord = nan(n,size(x,2));
for i = 1:n
    V1 = cord1(i,:);
    V2 = cord2(i,:);
    a = V2 - V1;
    b = x(i,:) - V1;
    aSquared = sqrt(dot(a,a));
    bSquared = sqrt(dot(b,b));
    abSq = aSquared.*bSquared;
    t  = dot(a,b)./(abSq);
%     ii = t>=0 & t<=1; %Angle between 0 and 90
%     theta = acos(t);
    bb = t.*bSquared;
    relDis(i,1) = bb./aSquared;
    projectedCoord(i,:) = a + bb;
end


end

