% Function to find the lowest error

function [erro,indx]=lowerror(dm,k)

no_feature=size(dm,1);
HIGH=9999;

for i=1:no_feature,
   if dm(i,i)==1
      kd(i)=HIGH;
   else
      
   if i==1
      dd=[HIGH,dm(1,2:no_feature)];
   elseif i==no_feature
      dd=[dm(1:no_feature-1,no_feature)',HIGH];
   else
      dd=[dm(1:i-1,i)',HIGH,dm(i,i+1:no_feature)];
   end
   for l=1:no_feature,
      if dm(l,l)==1
         dd(l)=HIGH;
      end
   end
   
   dd=sort(dd);
   kd(i)=dd(k);
   end

end


[erro,indx]=min(kd);
