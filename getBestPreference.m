function [ yp, xp,h ] = getBestPreference( x,y,pl )
% %getBestPreference Perform knee point detection
% to get the best clusters
% "Knee Point Detection in BIC for Detecting the Number of Clusters"
% Input:
%       x - X axis values
%       y - y axis values
%       pl - toggle plotting option
if(nargin==2)
    pl = false;
end
yp = 1;
xp = find(y==1);
% xp = xp(end);
if(size(x,1) ~= size(y,1))
    disp('ERROR')
    return;
end
ys = y;
a = 1;
b = .25*ones(4,1);
y = filter(a,b,y);
% y = smooth(y,'loess');
      
maxabd = abs(y(3:end)+y(1:end-2)-2.*y(2:end-1));
% Sort maxabd with values and indices
ix = zeros(size(maxabd,1),1);
uMaxabd = sort(unique(maxabd),'descend');
cnt = 0;
for i = 1:numel(uMaxabd)
    ii = find(maxabd==uMaxabd(i));
    ix(cnt+1:cnt+numel(ii)) = sort(ii,'descend');
    cnt = cnt+numel(ii);
end

n = 5;
ix = ix(1:n);
mangle = zeros(n,1);
for i = 1:n
    if(ix(i)>1)
        mangle(i) = atan(1/abs(y(ix(i)-1)-y(ix(i)))) + atan(1/abs(y(ix(i)+1)-y(ix(i))));   
    end
end
maxMangle = min(mangle);
uI = mangle==maxMangle;
% [mangle,i] = sort(mangle,'descend');
% uI = unique()
im = max(ix(uI));
y = ys;
yp = y(im-1);
xp = x(im-1);
h=[];

% plot figure
if(pl)
    h=figure;hold on;
    plot(x,y,'-r','MarkerFaceColor','r','MarkerEdgeColor','none');
    plot(xp,yp,'-ob','MarkerFaceColor','b','MarkerEdgeColor','none');
    hold off;
    ylabel('Number of clusters');xlabel('Preference')
    legend({'#Clusters';'Optimal Cluster'});
    xCent = min(x) + (max(x) - min(x))/2;
    yCent = (y(end) - y(1))/2;
    text(xCent,yCent,['Estimated Optimal Cluster -- ' num2str(yp)],'Color',[0 0 1],...
        'Fontweight','bold');
end

end

