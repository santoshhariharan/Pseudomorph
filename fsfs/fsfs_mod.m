% Reduces a feature set using the method described in:
% P. Mitra, C. A. Murthy and S. K. Pal, Unsupervised Feature Selection using 
% Feature Similarity, IEEE Transactions on Pattern Analysis and Machine
% Intelligence, Vol .24, No. 4, pp 301-312, April 2002. 

% Inputs:
% data : Data matrix containing the original feature set. Each column
% represents a feature.
% original_size : Number of features in the original data set.
% k : Scale parameter which decides the size of the reduced feature set.
% Approximately, k = original_size - size of the reduced feature set.

% Outputs:
% redu : Reduced features. A vector containing the feature numbers of the 
% original feature set which are retained in the reduced set.
% fwt : feature weights of the features in redu. 

function [redu,fwt]=fsfs_mod(dm,original_size,k)

no_feature=original_size;
% no_data=size(data,1);

%fprintf(1,'No. of features = %d\n',no_feature);

kk=[];
while k > (no_feature-1),
fprintf(1,'Give a smaller value of k\n');
k=input('k=  ');
end

% method=3;
% 1 = Feature Similarity: Correlation Coeff 
% 2 = Feature Similarity: Linear Regression error
% 3 = Feature Similarity: Maximal Information Compression Index

% Form the inter-feature distance matrix (Upper Triangular)
% fprintf(1,'Computing Feature Similarities..\n');
% for i=1:no_feature,
%    fprintf(1,'Similarity Computed for Feature No. %d\n',i);
%    for j=1:no_feature,
%       x1=data(:,i);x2=data(:,j);
%       if i < j
%          dm(i,j)=f2f(x1,x2,method);
%       else
%          dm(i,j)=0.0;
%       end
%    end
% end

drift=1.0;


% Form a vector containing the distance of k-NN for each feature.
for i=1:no_feature,
   if i==1
      dd=dm(1,2:no_feature);
   elseif i==no_feature
      dd=dm(1:no_feature-1,no_feature)';
   else
      dd=[dm(i,i+1:no_feature),dm(1:i-1,i)'];
   end
   dd=sort(dd);
   kd(i)=dd(k);
end
kd0=kd; % Store the original r_k's


% Condense the feature set
fprintf(1,'\nClustering the Features..\n');
rfs=[];rfd=[];ee=[];dmt=dm;lower=9999;iter=0;prev_lower=9999;
tagm=ones(1,no_feature);
while (no_feature-trace(dm)) >0,
   iter=iter+1;
   if k > (no_feature-trace(dm)-1)
      k=no_feature-trace(dm)-1;
   end
   if k<=0
      break;
   end
   prev_lower=lower;
   [lower,fetr]=lowerror(dm,k);
  
   
   % Adjust k
   while lower > drift*prev_lower,
      k=k-1;
      if k==0
         break;
      end
      [lower,fetr]=lowerror(dm,k);
   end
   
   if k <=0
      break;
   end
   dm=updatedm(dm,fetr,k);
          
   kk=[kk;k];
   ee=[ee;lower];
     
   tagm(fetr)=0;
   for i=1:no_feature,
      for j=1:no_feature,
         if dm(i,i)==1
            tagm(i)=0;
         end
      end
   end
   
end

for i=1:no_feature,
   if dm(i,i)==0
      rfs=[rfs;i];
      rfd=[rfd;kd0(i)];
   end
end

fprintf(1,'Features Clustered.\n');
     
redu=rfs;
fwt=rfd;








