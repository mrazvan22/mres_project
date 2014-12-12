%% Setting up simulated data
clear all
close all

nr_roi = 4;
nr_pat = 100;

order_events_sim1 = 1:nr_roi;
order_events_sim2 = [1 3 2 4];
atrophy_model_sim1 = zeros(nr_roi, nr_roi);
atrophy_model_sim2 = zeros(nr_roi, nr_roi);
for roi = 1:nr_roi,
    
    atrophy_model_sim1(order_events_sim1==roi, roi:end) = 1;
    atrophy_model_sim2(order_events_sim2==roi, roi:end) = 1;
    
end

I_rand = ceil(nr_roi*rand(nr_pat, 1));
p_atrophy = zeros(nr_roi, nr_pat);
for pat = 1:nr_pat,
    
    if rand < 0.5,
        
        p_atrophy(:, pat) = atrophy_model_sim1(:, I_rand(pat));
        
    else
        
        p_atrophy(:, pat) = atrophy_model_sim2(:, I_rand(pat));
        
    end
    
end
data_atrophy = p_atrophy;
I_randperm = randperm(nr_roi*nr_pat);
data_atrophy(I_randperm(1:10)) = 1-data_atrophy(I_randperm(1:10));

%% Performing mcmc
c_lim = [0.001 0.1];
d_lim = [0.1 0.8];
nr_it_mcmc = 1e4;
nr_it_burnin = 1e4;
nr_it_hillclimb = 1e3;
nr_hillclimb = 5;
thinning = 1;

[parm_struct, diag_struct] = ...
    EventOnsetModel(data_atrophy, nr_it_hillclimb, ...
    nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb);

hist2_mat = zeros([nr_roi nr_roi]);
figure(1), clf

thr_vec = [1e-4 1e-2 0.1 0.2 0.4 0.6];
for thr = 1:length(thr_vec),
    
    p_A_D = data_atrophy;
    p_A_D(p_A_D < thr_vec(thr)) = thr_vec(thr);
    p_A_D(p_A_D > 1-thr_vec(thr)) = 1-thr_vec(thr);
        
    [parm_struct{thr}, diag_struct{thr}] = ...
        AtrophyModelMCMC2a(p_A_D, ...
        nr_it_hillclimb, nr_it_burnin, nr_it_mcmc, ...
        thinning, nr_hillclimb, 2);
    
end

[ord inds_av] = sort(mean(parm_struct.order_events, 2));
for roi1 = 1:nr_roi,
    
    for roi2 = 1:nr_roi,
        
        hist2_mat(roi1, roi2) = sum(parm_struct.order_events(inds_av(roi1), :) == roi2);
        
    end
    
end
imagesc(hist2_mat)

eval(sprintf('save simRegions%dPatients%d', nr_roi, nr_pat))

