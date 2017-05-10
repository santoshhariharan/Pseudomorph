function C = clsIn(data,beta,~)



if(isempty(data))
    disp('Datais empty');
    return;
end
% Initialize 
C.minClsSize = 5;
C.maxCls = 10;
C.minCls = 1;
C.S = [];
C.pmin = 0;
C.pmax = 0;
C.pmed = 0;

if(nargin == 1)
    beta = .05;
%     type=1;
% elseif(nargin==2)
%     if(beta <=0 || beta >=1)
%         beta = .05;
%     end
%     type = 1;
% else
%     if(strcmpi(type,'Cosine'))
%         type = 3;
%     elseif(strcmpi(type,'Gaussian'))
%         type = 2;
%     else
%         type = 1;
%     end
end
% type=2;
% [mm] = size(data,1);
% M=mm*mm-mm;
% sim=zeros(M,3); % Make ALL N^2-N similarities
% j=1;
% if(type==1 || type ==2)
% for i=1:mm
%   for k=[1:i-1,i+1:mm]
%     sim(j,1)=i; sim(j,2)=k; sim(j,3)=-sum((data(i,:)-data(k,:)).^2);    
%     j=j+1;
%   end;
% end;
% if(type==2)
%     sim(:,3) = (1/sqrt(2*pi*2))*exp(sim(:,3)./8);
% end
% else
%     for i=1:mm
%       for k=[1:i-1,i+1:mm]
%           sim(j,1)=i; sim(j,2)=k;
%           sim(j,3) = sum(data(i,:).*data(k,:));
%           tmp = sqrt(sum(data(i,:).^2)).*sqrt(sum(data(k,:).^2));
%           sim(j,3) = sim(j,3)./tmp;
%           j=j+1;
%       end
%     end
%     
% end
sim = -1*pdist2(data,data,'Euclidean');
% sim = -1*sqDistance(data',data');
x_x = logical(tril(ones(size(sim,1),size(sim,1)),-1));
% assignin('base','x_x',x_x)
% sim = -pdist2(data,data,'Euclidean');

C.pmed = median(sim(x_x)); 
[C.pmin, C.pmax] = preferenceRange(sim);
% C.minClsSize = round(beta*mm);
% C.maxCls = round(mm./C.minClsSize);
C.S = sim;
end