
a1 = [2 3;3 5; 1 6;4 3];
% plot(a1(:,1),a1(:,2),'*r');

a2 = [0 0 1 1;0 0 1 1; 0 1 0 0;1 1 0 0];
d = a2.*pdist2(a1,a1);
[a,b] = computePathlength(a2,w);
