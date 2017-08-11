function [ D,N ] = computePathlength( A, W )
%computePathlength Computes path length of all points in the graph to all
%                   other points using Dijkstra's algorithm

% A - Adjacency matrix - m x m where m is the number of vertices
% W - Weights for each of the edge
% D - Path length - m x m matrix
% N - Number of nodes traversed along the path


m = size(A,1);
[p,q] = size(W);
D = inf(m,m);
N = inf(m,m);
if(p == q && p == m)
    w = W.*A;
else
    w = zeros(m,m);
    for i = 1:p
        w(W(i,1),W(i,2)) = W(i,3);
    end
    w = w + w';
end
jj = w==0;
w(jj) = inf;
fprintf('Computing path length...............')
for iVertex = 1:m % For each vertex
    fprintf('\b\b\b\b\b\b\b\b%7.3f%%',iVertex*100/m);
    % Initialize
%     disV = inf(1,m);
    prevV = zeros(m,1);
    
    Q = true(1,m);
    Q(iVertex) = false;    
    disV = w(iVertex,:);disV(iVertex) = 0;
%     prevV = zeros(m,1);    
%     allW = w(iVertex,:);
    while sum(Q)>0
        tmp = find(Q);
        [du,u] = min(disV(1,Q));
        u = tmp(u); 
        Q(u) = false; % Remove u from Q
        tmpX = A(u,:);
        allU = w(u,:); % Find neighbors of u        
        ii = find(logical(tmpX)); % Neighbour of u Connected edge that is not visited        
        if(isempty(ii))
%             allU = w(:,u)';
%             ii = find(allU ~= inf & Q);
%             if(isempty(ii))
                disp('breaking');
                break;
%             end
        end
        for j = 1:numel(ii)
            alt = du + allU(ii(j));
            if(alt<disV(ii(j)))
                disV(ii(j)) = alt;
                prevV(ii(j)) =u;
            end            
        end
    end
    D(iVertex,:) = disV;
    N(iVertex,:) = prevV';
end

% Reset directly connected euclidean distances
% for i = 1:p
%     D(W(i,1),W(i,2)) = W(i,3);
% end



end
