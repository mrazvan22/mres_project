function logLik = logLikelihoodAtrophyModelPuolamaki(data_atrophy, ...
    atrophy_model, c, d)

[nr_roi, nr_pat] = size(data_atrophy);
logLik = zeros(nr_pat, 1);
for pat = 1:nr_pat,
    
    p_alpha = zeros(nr_roi, 1);
    p_beta = zeros(nr_roi, 1);
    p_gamma = zeros(nr_roi, 1);
    p_delta = zeros(nr_roi, 1);
    for roi = 1:nr_roi,
                        
        I_alpha = find((data_atrophy(:, pat) == 1) & ...
            (atrophy_model(:, roi) == 0));
        p_alpha(roi) = length(I_alpha)*log(c);
        I_beta = find((data_atrophy(:, pat) == 1) & ...
            (atrophy_model(:, roi) == 1));
        p_beta(roi) = length(I_beta)*log(1-c);
        I_gamma = find((data_atrophy(:, pat) == 0) & ...
            (atrophy_model(:, roi) == 1));
        p_gamma(roi) = length(I_gamma)*log(d);
        I_delta = find((data_atrophy(:, pat) == 1) & ...
            (atrophy_model(:, roi) == 1));
        p_delta(roi) = length(I_delta)*log(1-d);
        %         I_alpha = find((data_atrophy(:, pat) == 1) & ...
        %             (atrophy_model(:, roi) == 0));
        %         p_alpha(roi) = sum(log(p_noatrophy(I_alpha, pat)));
        %         I_beta = find((data_atrophy(:, pat) == 1) & ...
        %             (atrophy_model(:, roi) == 1));
        %         p_beta(roi) = sum(log(1 - p_noatrophy(I_beta, pat)));
        %         I_gamma = find((data_atrophy(:, pat) == 0) & ...
        %             (atrophy_model(:, roi) == 1));
        %         p_gamma(roi) = sum(log(1 - p_noatrophy(I_gamma, pat)));
        %         I_delta = find((data_atrophy(:, pat) == 0) & ...
        %             (atrophy_model(:, roi) == 0));
        %         p_delta(roi) = sum(log(p_noatrophy(I_delta, pat)));
        
    end
    logLik_pat = p_alpha + p_beta + p_gamma + p_delta;
    I_max = find(logLik_pat == max(logLik_pat));
    logLik(pat) = logLik_pat(I_max(1));
    
end
logLik = sum(logLik);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------