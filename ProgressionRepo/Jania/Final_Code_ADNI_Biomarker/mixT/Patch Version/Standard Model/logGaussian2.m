function logProb= logGaussian2(data,mean,covar)

nDim = size(data,1);
covar(covar==0)=0.1;

if (size(covar,2)==1)%if cov is diagonal
    logProb = -0.5*nDim*log(2*pi)-0.5*sum(covar)-0.5*sum(((data-mean).^2)./covar,1);
else %if it is general (slower)
    logProb= -0.5*nDim*log(2*pi)-0.5*log(det(covar))-0.5*(data-mean)'*inv(covar)*(data-mean);
end
