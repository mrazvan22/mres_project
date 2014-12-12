function rt=mixTDemo(data,dataIDMatrix,N_T_EST)
%demonstration program that fits N gaussians in M dimensional space

close all;





%cov type parameter
%covType = 0;  %unconstrained - individual covariance matrices
%covType = 1;  %unconstrained - shared covariance matrices
%covType = 2;  %diagonal - individual covariance matrices
covType = 3;  %diagonal - shared covariance matrices
%covType = 4;  %spherical - individual covariance matrices
%covType = 5;  %spherical - shared covariance matrices
 
%tEst = mixTFit(data,N_T_EST,covType);
tEst = mixTFitMosaicModelFast(data,N_T_EST,covType,dataIDMatrix);

% 
% fprintf('Estimated Gauss 1\n');
% m = tEst(1).mean
% c = tEst(1).cov
% p = tEst(1).prior
% d = tEst(1).dof
% 
% fprintf('Estimated Gauss 2\n');
% 
% tEst(2).mean
% tEst(2).cov
% tEst(2).prior
% tEst(2).dof

 
rt=tEst;