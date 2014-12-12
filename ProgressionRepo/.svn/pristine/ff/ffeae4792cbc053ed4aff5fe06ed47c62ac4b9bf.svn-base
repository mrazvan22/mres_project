randn('state', 0); rand('state', 0);
ndim = 3;
gmix = gmm(ndim, 2, 'spherical');
ndat1 = 300; ndat2 = 300; ndata = ndat1+ndat2;
gmix.centres =  [0.75 0.75 0.75
    0.25 0.25 0.25];
gmix.covars = repmat(0.01, 1, ndim);
x = gmmsamp(gmix, ndata);

for ncentres = 1:10,
    
    opt_mix = foptions;
    opt_mix(14) = 1e4;
    mix = gmm(ndim, ncentres, 'spherical');
    mix = gmminit(mix, x, opt_mix);
    [mix, opt_mix] = gmmem(mix, x, opt_mix);
    
    logLik = -opt_mix(8);
    npar = ncentres - 1 + ncentres*(ndim + 1);
    aic(ncentres) = aic_hf(npar, ndata, logLik);
    
end
plot(aic)