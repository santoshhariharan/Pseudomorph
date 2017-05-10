function [ opFeatures,featureScores ] = unsupervisedPCASelect( X,numfeatures )
%unsupervisedPCASelect Summary of this function goes here
%   Detailed explanation goes here

if(numfeatures>size(X,2))
    opFeatures = true(1,size(X,2));
    return;
end
[COEFF,~,latent] = princomp(X);

latentCumSum = cumsum(latent./sum(latent));

numPC = sum(latentCumSum<=.7);
if(numPC>15)
    numPC = 15;
end
COEFF = COEFF(:,1:numPC);
latent = latent(1:numPC);
featureScores = sum(bsxfun(@times,COEFF,latent'),2);
[~,I] = sort(featureScores,'descend');
I = I(1:numfeatures);
opFeatures = false(1,size(X,2))';
opFeatures(I) = true;

end

