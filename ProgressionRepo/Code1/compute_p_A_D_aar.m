function p_A_D = compute_p_A_D_aar(atrophy_patient, atrophy_control, ...
    atrophy_rate_patient, atrophy_rate_control, diff_t_patient, diff_t_control)

[nr_roi, nr_pat] = size(atrophy_patient);
nr_con = size(atrophy_control, 2);

X_control = [ones(nr_con, 1) diff_t_control];
X_patient = [ones(nr_pat, 1) diff_t_patient];
p_noA_D1 = zeros(nr_roi, nr_pat);
for roi = 1:nr_roi,
    
    b = X_control\atrophy_control(roi, :)';
    e = atrophy_control(roi, :)' - X_control*b;
    sig = std(e);
    p = normpdf(atrophy_patient(roi, :)', X_patient*b, sig);
    p(atrophy_patient(roi, :)' > X_patient*b) = 1-eps;
    p_noA_D1(roi, :) = p;
    
end

p_noA_D2 = zeros(nr_roi, nr_pat);
for roi = 1:nr_roi,
    
    [mu, sig] = normfit(atrophy_rate_control(roi, :));
    p = normpdf(atrophy_rate_patient(roi, :)', mu, sig);
    p(atrophy_rate_control(roi, :) > mu) = 1-eps;
    p_noA_D2(roi, :) = p;
    
end
p_noA_D = min(cat(3, p_noA_D1, p_noA_D2), [], 3);
p_A_D = eps*ones(nr_roi, nr_pat);
I_sig = find(p_noA_D < 0.05);
p_A_D(I_sig) = 1 - p_noA_D(I_sig);

% % -------------------------------------------------------------------------
% % -------------------------------------------------------------------------
% function [offset, r, sig] = mcmc_regress(X, y)
% 
% nr_it_burnin = 1e4;
% nr_it_mcmc = 1e4;
% X_reg = [ones(length(X), 1) X];
% b = X_reg\y;
% if (b(2)) > 0,
%     
%     r_current = -b(2);
%     
% else
%     
%     r_current = b(2);
%     
% end
% sig_current = std(y - X_reg*b);
% offset_current = b(1);
% logLik_current = sum(log(normpdf(y, offset_current - r_current*X, sig_current)));
% std_update_r = abs(r_current)*0.05;
% std_update_sig = sig_current*0.1;
% std_update_offset = abs(offset_current)*0.05;
% 
% offset = zeros(nr_it_mcmc, 1);
% r = zeros(nr_it_mcmc, 1);
% sig = zeros(nr_it_mcmc, 1);
% cnt_it = 1:500:(nr_it_burnin + nr_it_mcmc);
% for it = 1:(nr_it_burnin + nr_it_mcmc),
%     
%     flg_continue =1;
%     while flg_continue,
%         
%         r_new = normrnd(r_current, std_update_r);
%         if r_new < 0
%             
%             flg_continue = 0;
%             
%         end
%         
%     end
%     flg_continue = 1;
%     while flg_continue,
%         
%         sig_new = normrnd(sig_current, std_update_sig);
%         if sig_new > 0,
%             
%             flg_continue = 0;
%             
%         end
%         
%     end
%     offset_new = normrnd(offset_current, std_update_offset);
%     Lik_new = sum(log(normpdf(y, offset_new - r_new*X, sig_new)));
%     alpha = Lik_new/Lik_current;
%     if alpha > rand,
%         
%         r_current = r_new;
%         sig_current = sig_new;
%         offset_current = offset_new;
%         
%     end
%     if it > nr_it_burnin,
%         
%         r(it-nr_it_burnin) = r_current;
%         sig(it-nr_it_burnin) = sig_current;
%         offset(it-nr_it_burnin) = offset_current;
%         
%     end
%     
% end
% 
% 
% 
