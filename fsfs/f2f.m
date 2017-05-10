% Function to compute correlation between two variables

function dist=f2f(x1,x2,method)

no_x1=size(x1,1);
no_x2=size(x2,1);
dist=0.0;

% Distance = Correlation Cefficient
if method==1
   num=0.0;den1=0.0;den2=0.0;
   x1bar=mean(x1);x2bar=mean(x2);
   for i=1:no_x1,
      num=num+abs(x1(i)*x2(i)-x1bar*x2bar);
      den1=den1+(x1(i)-x1bar)^2;den2=den2+(x2(i)-x2bar)^2;
   end
   dist=num/sqrt(den1*den2);
  
      
   % Distance = regression error
elseif method==2
   num=0.0;den=0.0;x1bar=mean(x1);x2bar=mean(x2);   
   for i=1:no_x1,
      num=num+x1(i)*x2(i)-x1bar*x2bar;
      den=den+x1(i)^2-x1bar^2;
   end
   a=num/den;b=x2bar-a*x1bar;
   
   dist=0.0;
   
   for i=1:no_x1,
      dist=dist+abs(x2(i)-a*x1(i)-b)/sqrt(a^2+b^2);
   end
   dist=dist/no_x1;
   
elseif method==3
   sxy=0.0;sx=0.0;sy=0.0;mnx1=0.0;mnx2=0.0;
   for i=1:no_x1,
      sxy=sxy+x1(i)*x2(i);
      sx=sx+x1(i)^2;
      sy=sy+x2(i)^2;
      mnx1=mnx1+x1(i);
      mnx2=mnx2+x2(i);
   end
   mnx1=mnx1/no_x1;
   mnx2=mnx2/no_x2;
   
   sxy=(sxy/no_x1)- mnx1*mnx2;
   sx=(sx/no_x1)-mnx1^2;
   sy=(sy/no_x1)-mnx2^2;
   
     if (sx-sy) ==0
      theta=0.5*pi/2;
      else
         theta=0.5*atan(2*sxy/(sx-sy));
      end
      
   a=-cot(theta);
   b=mean(x1)*cot(theta)+mean(x2);
   
   dist=0.0;
   for i=1:no_x1,
      dist=dist+abs(x2(i)-a*x1(i)-b)/sqrt(a^2+b^2);
   end
   dist=dist/no_x1;  
   
elseif method==4
   cv=cov(x1,x2)/no_x1;
   dist=min(eig(cv));
   
   
   end
