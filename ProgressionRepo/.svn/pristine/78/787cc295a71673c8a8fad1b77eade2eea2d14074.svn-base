function rt=version3FindSingleTermsTDist(X,library)

N_IND=size(library,1);
N_APP=size(library,2);


%sum over H Model1

for(cInd=1:N_IND)
    %sum over A Model1
    for(cApp=1:N_APP)
        prAM1(cApp)=getLogTProb(X,library(cInd,cApp).mean,library(cInd,cApp).cov,library(cInd,cApp).dof);
    end
    prH1M1(cInd)=sumLog(prAM1);
end

rt=sumLog(prH1M1);


%========================================================================
function rt=sumLog(lst)
% a function that takes a list of terms, rescales the logs and return the
dim=size(lst,2);
logSumTerms=0;
sumTerms=0;
%find the rescale factor
rescaleFactor=max(lst(:));
%subtract the rescale factor
scaledTerm=lst-rescaleFactor;
%takes the exp
scaledTerm=exp(scaledTerm);
%add them together
sumTerms=sum(scaledTerm(:));
%convert back to log  probabilities and add the rescale factor
logSumTerms=log(sumTerms)+rescaleFactor;
 
rt=logSumTerms;

%========================================================================
function logProb= getLogMVarGaussProb(data,mean,covar)

nDim = size(data,1);
covar(covar==0)=realmin;

if (size(covar,2)==1)%if cov is diagonal
    logProb = -0.5*nDim*log(2*pi)-0.5*sum(covar)-0.5*sum(((data-mean).^2)./covar,1);
else %if it is general (slower)
    logProb= -0.5*nDim*log(2*pi)-0.5*log(det(covar))-0.5*(data-mean)'*inv(covar)*(data-mean);
end;