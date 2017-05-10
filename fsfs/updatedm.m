% Function to recompute the distance matrix during clustering

function dm1=updatedm(dm,indx,k)
no_feature=size(dm,1);

HIGH=9999;
i=indx;

   if i==1
      dd=[HIGH,dm(1,2:no_feature)];
      for l=1:no_feature,
         if dm(l,l)==1
            dd(l)=HIGH;
         end
      end
      
   elseif i==no_feature
      dd=[dm(1:no_feature-1,no_feature)',HIGH];
      for l=1:no_feature,
         if dm(l,l)==1
            dd(l)=HIGH;
         end
      end

   else
      dd=[dm(1:i-1,i)',HIGH,dm(i,i+1:no_feature)];
      for l=1:no_feature,
         if dm(l,l)==1
            dd(l)=HIGH;
         end
      end

   end
      
   [dd,dindx]=sort(dd);
   
   dm1=dm;
      
   for l=1:k,
      indx1=dindx(l);
      dm1(indx1,indx1)=1;
   end
   
      
   
   


   
   
