
a1 = [2 3;3 5; 1 6;4 3];
% plot(a1(:,1),a1(:,2),'*r');

a2 = [0 0 1 1;0 0 1 1; 0 1 0 0;1 1 0 0];
a2 = a2+a2';
a2 = double(a2>0);
d = a2.*pdist2(a1,a1);
[a,b] = computePathlength(a2,d);
cc = dijkstra(a2,d);
%% 
    n = 10; A = zeros(n); xy = 10*rand(n,2);
    tri = delaunay(xy(:,1),xy(:,2));
    I = tri(:); J = tri(:,[2 3 1]); J = J(:);
    IJ = I + n*(J-1); A(IJ) = 1;
    a = (1:n); b = a(ones(n,1),:);
    C = round(reshape(sqrt(sum((xy(b,:) - xy(b',:)).^2,2)),n,n));
    %%
    A = A+A';
    A = double(A>0);
    C = C+C';
    C = C/2;
    [costs,paths] = dijkstra(A,C);
    [a,b] = computePathlength(A,C);
    isequal(a,costs)
%     
