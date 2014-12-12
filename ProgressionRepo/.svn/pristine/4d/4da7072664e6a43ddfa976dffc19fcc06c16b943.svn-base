function r= mixTGenerate(tDist,nData)
%generate data from mixture of Gaussians

N_DIST = length(tDist);
N_DIM = size(tDist(1).mean,1);

%generate data from gaussians
prTDistCum = cumsum([tDist(:).prior]);

%decide if covariance is diagonal or not
if (size(tDist(1).cov,2)>1)
    DIAG_FLAG = 0;
else
    DIAG_FLAG = 1;
end;

%do cholesky decomp of cov matrix
if (DIAG_FLAG)
    for (cDist = 1:N_DIST)
        tDist(cDist).stdDev = sqrt(tDist(cDist).cov);     
    end;
else
    for (cDist = 1:N_DIST)
         tDist(cDist).cholCov = chol(tDist(cDist).cov);
    end;
end;


r = zeros(N_DIM,nData);

for (cData = 1:nData)
    %choose which Gaussian
    randVal = rand(1);
    decision = prTDistCum-randVal;
    %retain gaussians where cum prob is greater
    decision(find(decision<0))=1.1;
    %choose closest Gaussian
    distChosen = find(decision == min(decision));
    %generate this point and store
    uVal = gamrnd(tDist(distChosen).dof/2,tDist(distChosen).dof/2,1,1);
    covScaling = 1/sqrt(uVal);
    if (DIAG_FLAG)
        r(:,cData) = tDist(distChosen).mean+covScaling*tDist(distChosen).stdDev.*randn(N_DIM,1);        
    else
        r(:,cData) = tDist(distChosen).mean+covScaling*tDist(distChosen).cholCov*randn(N_DIM,1);
    end
end;
