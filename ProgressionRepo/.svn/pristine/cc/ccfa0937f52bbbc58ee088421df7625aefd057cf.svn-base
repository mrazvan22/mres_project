function rt=mixGaussPatchBasedVersion3(data,dataIDMatrix,N_IND_EST,N_APP_EST)
%demonstration program that fits N gaussians in M dimensional space

close all;

data=squeeze(data)';


%cov type parameter
%covType = 0;  %unconstrained - individual covariance matrices
%covType = 1;  %unconstrained - shared covariance matrices
%covType = 2;  %diagonal - individual covariance matrices
covType = 3;  %diagonal - shared covariance matrices
%covType = 4;  %spherical - individual covariance matrices
%covType = 5;  %spherical - shared covariance matrix


% data=data(:,1:8);
% dataIDMatrix=dataIDMatrix(1:8,1:2);


N_ATTEMPTS = 1;
%tDist= mixTFitAppearanceLog(data,N_IND_EST,N_APP_EST,covType,N_ATTEMPTS,dataIDMatrix);
tDist= mixTFitAppearanceModelFast(data,N_IND_EST,N_APP_EST,covType,N_ATTEMPTS,dataIDMatrix);

%save([destinationDir 'version3Gaussian32i2a.mat'],'gaussEst','prData');


rt=tDist;
