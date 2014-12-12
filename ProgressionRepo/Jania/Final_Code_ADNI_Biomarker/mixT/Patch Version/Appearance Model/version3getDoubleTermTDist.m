function rt=version3getDoubleTerm(X, Xp,library)

N_IND=size(library,1);
N_APP=size(library,2);

% for(cInd=1:N_IND)
% 
%     for(cInd=1:N_IND)
%         %sum over A Model1
%         for(cApp=1:N_APP)
%             prAM1(cApp)=getLogMVarGaussProb(X,library(cInd,cApp).mean,library(cInd,cApp).cov)+...
%                 getLogMVarGaussProb(Xp,library(cInd,cApp).mean,library(cInd,cApp).cov);
%         end
%         prH1M1(cInd)=sumLog(prAM1);
%     end
% 
%     summedOverThisH=sumLog(prH1M1);
% 
%     sumH(cInd)=summedOverThisH;
% end
% 
% rt=sum(sumH);
% 

for(cInd=1:N_IND)
   
    %sum over A Model1
    for(cApp=1:N_APP)
        prXGivenHiAij(cInd,cApp)=getLogTProb(X,library(cInd,cApp).mean,library(cInd,cApp).cov,library(cInd,cApp).dof)+...
            log(library(cInd,cApp).prior);
    end
    prXGivenHi=sumLog(prXGivenHiAij(cInd,:));
    for(cApp=1:N_APP)
        prXPGivenHiAij(cInd,cApp)=getLogTProb(Xp,library(cInd,cApp).mean,library(cInd,cApp).cov,library(cInd,cApp).dof)+...
            log(library(cInd,cApp).prior);
    end

    prXPGivenHi=sumLog(prXPGivenHiAij(cInd,:));
    prH1M1(cInd)=prXGivenHi+prXPGivenHi;


end
     
%       %Marginalization alternative
       prHGivenData=sumLog(prH1M1);  
      
rt=prHGivenData;

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

% %========================================================================
% function logProb= getLogMVarGaussProb(data,mean,covar)
% 
% nDim = size(data,1);
% covar(covar==0)=realmin;
% if (size(covar,2)==1)%if cov is diagonal
%     logProb = -0.5*nDim*log(2*pi)-0.5*sum(covar)-0.5*sum(((data-mean).^2)./covar,1);
% else %if it is general (slower)
%     logProb= -0.5*nDim*log(2*pi)-0.5*log(det(covar))-0.5*(data-mean)'*inv(covar)*(data-mean);
% end;
