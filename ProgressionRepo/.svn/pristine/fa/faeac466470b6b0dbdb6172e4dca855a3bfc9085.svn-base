function logProb= getLogMVarGaussProb(data,mean,covar)
%return log probability of multivariate gaussian

[nDim nData] = size(data);
logProb = zeros(1,nData);

if (size(covar,2)==1)%if cov is diagonal
    for (cData = 1:nData)
        logProb(cData) = -0.5*nDim*log(2*pi)-0.5*sum(covar)-0.5*sum(((data(:,cData)-mean).^2)./covar,1);
    end;
else %if it is general (slower)
    invCovar = inv(covar);
    logDetCovar = log(det(covar));
    for (cData = 1:nData)
        logProb(cData)= -0.5*nDim*log(2*pi)-0.5*logDetCovar-0.5*(data(:,cData)-mean)'*invCovar*(data(:,cData)-mean);
    end
end;
