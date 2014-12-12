function logProb = mixGaussLogProb(data,mixGauss)

nGauss = length(mixGauss);
nData= size(data,2);

logProb = zeros(1,nData);

for (cGauss = 1:nGauss)
    logProb = logProb+repmat(log(mixGauss(cGauss).prior),1,nData)+...
        getLogMVarGaussProb(data,mixGauss(cGauss).mean,mixGauss(cGauss).cov);
end;            
 