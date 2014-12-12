nr_it_mcmc = 1e4;
nr_it_burnin = 1e4;
nr_it_hillclimb = 2e3;
thinning = 1e0;
nr_hillclimb = 10;

thr_vec = [1e-6 1e-1 0.2 0.4];
for it_thr1 = 1:length(thr_vec),
    
    for it_thr2 = 1:length(thr_vec),
        
        p_A_D_filt = p_A_D(:, :, 1);
        p_A_D_filt(p_A_D_filt < thr_vec(it_thr1)) = thr_vec(it_thr1);
        p_A_D_filt(p_A_D_filt > (1-thr_vec(it_thr2))) = 1-thr_vec(it_thr2);
        [parm_struct{it_thr1, it_thr2}, diag_struct{it_thr1, it_thr2}] = ...
            AtrophyModelMCMC2(p_A_D_filt, nr_it_hillclimb, ...
            nr_it_burnin, nr_it_mcmc, thinning, nr_hillclimb, version_mcmc);
        fprintf('it_thr1: %.3f\tit_trh2: %.3f\n', thr_vec(it_thr1), thr_vec(it_thr2))
        
    end
        
end
hist2_mat = zeros([nr_events nr_events length(thr_vec) length(thr_vec)]);
for it_thr1 = 1:length(thr_vec),
    
    for it_thr2 = 1:length(thr_vec),
        
        [ord inds_av] = sort(mean(parm_struct{it_thr1, it_thr2}.order_events, 2));
        for roi1 = 1:nr_events,
            
            for roi2 = 1:nr_events,
        
                hist2_mat(roi1,roi2, it_thr2, it_thr1) = ...
                    sum(parm_struct{it_thr1, it_thr2}.order_events(inds_av(roi1), :) == roi2);
                
            end
            
        end
        
    end
    
end

for it_thr1 = 1:length(thr_vec),
    
    for it_thr2 = 1:length(thr_vec),
        
        logLik_max(it_thr1, it_thr2) = diag_struct{it_thr1, it_thr2}.logp_order_events_max;
        
    end
end
figure(1), clf
imagesc(logLik_max)

cnt_fig = 1;
figure(2), clf
for it_thr1 = 1:length(thr_vec),
    
    for it_thr2 = 1:length(thr_vec),
        
        subplot(3, 3, cnt_fig)
        imagesc(hist2_mat(:, :, it_thr2, it_thr1))
        title(sprintf('thr1: %d\tthr2: %d', it_thr1, it_thr2))
        cnt_fig = cnt_fig + 1;
        
    end
    
end
        
        
        
        
