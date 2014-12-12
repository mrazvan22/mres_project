 
function logProb= getLogTProb(data,meanVal,covar,dof)

 [N_DIM N_DATA] = size(data);
 xMinusMean =data-repmat(meanVal,1,N_DATA);
 
 
 if (size(covar,2)==1)%if cov is diagonal
     deltaVal = sum(xMinusMean.*(repmat(1./covar,1,N_DATA)).*xMinusMean,1);
     detCov = prod(covar);
 else %if it is general (slower)
     deltaVal = sum((xMinusMean'*inv(covar))'.*xMinusMean,1);
     detCov = det(covar);
 end;
 
 if(detCov==inf)
     detCov=realmax;
 end

 
 if(detCov==0)
     detCov=realmin;
 end
 
 logProb = log(gamma(dof/2+N_DIM/2))-log(gamma(dof/2));
 logProb = logProb-0.5*log(detCov);
 logProb = logProb-(N_DIM/2)*log(dof*pi);
 logProb = logProb-(dof/2+N_DIM/2)*log(1+deltaVal/dof);