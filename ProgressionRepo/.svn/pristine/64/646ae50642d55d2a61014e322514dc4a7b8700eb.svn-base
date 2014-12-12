clear all
close all

data = [mvnormrnd(repmat([-5 5 0]', 20, 1), 1, 4000);...
    mvnormrnd(repmat([5 0 -5]', 20, 1), 0.5, 4000)];

figure(1),
scatter3(data(:, 1), data(:, 2), data(:, 3))

cfg = [];
cfg.gmmoutput = 'best';
cfg.gmmsearch = 'bic';
cfg.minkernel = 1;
cfg.maxkernel = 20;
cfg.maxemiter = 200;

gmm = gmmfit(data, cfg);
    