function mixTDemo
%demonstration program that fits N gaussians in M dimensional space

close all;


%create data
N_DIST = 2;
N_DIM = 1;

%set up prior probability of t-distributions
prDist = 0.2+0.8*rand(N_DIST,1);
prDist = prDist/sum(prDist);

%fill in actual mean, prior and covariance
for (c1 = 1:N_DIST)
    tDist(c1).mean = 17*rand(N_DIM,1);
    tDist(c1).prior = prDist(c1);
    tDist(c1).dof = rand(1)*5+1.0;
    nDimNoise = randperm(N_DIM);
    F = randn(N_DIM,nDimNoise(1));
    tDist(c1).cov = F*F'+0.1*eye(N_DIM)+diag(0.2*rand(N_DIM,1));
end;


%generate the data
N_DATA = 1000;
  data = mixTGenerate(tDist,N_DATA);
%  save('data.mat','data');
%load('data.mat','data');

%fit gaussians to data
N_T_EST = 2;

%cov type parameter
%covType = 0;  %unconstrained - individual covariance matrices
%covType = 1;  %unconstrained - shared covariance matrices
%covType = 2;  %diagonal - individual covariance matrices
covType = 3;  %diagonal - shared covariance matrices
%covType = 4;  %spherical - individual covariance matrices
%covType = 5;  %spherical - shared covariance matrices

dataIDMatrix=createConseqDataID(2, 1000, 500);
tEst = mixTFit(data,N_T_EST,covType);
%tEst = mixTFitMosaicModelFast(data,N_T_EST,covType,dataIDMatrix);


fprintf('Estimated Gauss 1\n');
m = tEst(1).mean
c = tEst(1).cov
p = tEst(1).prior
d = tEst(1).dof

fprintf('Estimated Gauss 2\n');

tEst(2).mean
tEst(2).cov
tEst(2).prior
tEst(2).dof

fprintf('Actual Gauss 1\n');

tDist(1).mean
tDist(1).cov
tDist(1).prior
tDist(1).dof

fprintf('Actual Gauss 2\n');

tDist(2).mean
tDist(2).cov
tDist(2).prior
tDist(2).dof

i
