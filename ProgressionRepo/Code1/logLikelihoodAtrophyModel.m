% [logLik, logLikmat, model_bestfit] = logLikelihoodAtrophyModel(p_A_D, ...
%     order_events, version)
% Version signifies the way likelihood is calculated
% version == 1: the stage with the highest likelihood for each subject is
% taken
% version == 2: The likelihood of all stages is summed for each subject,
% then the product over all subjects is taken
% version == 3 The likelihood for all patients is summed for each stage,
% then the product over all stages is taken

function [logLik, logLikmat, model_bestfit] = logLikelihoodAtrophyModel(p_A_D, ...
    order_events, version)

[nr_roi, nr_pat] = size(p_A_D);
atrophy_model = zeros(nr_roi, nr_roi);
for roi = 1:nr_roi,
    
    atrophy_model(order_events==roi, roi:end) = 1;
    
end

logp_A_D = log(p_A_D);
p_noA_D = 1- p_A_D;
p_noA_D(p_noA_D == 0) = eps;
logp_noA_D = log(p_noA_D);
atrophy_model = atrophy_model == 1;
noatrophy_model = atrophy_model == 0;

if ~exist('version', 'var'),
    
    version = 2;
    
end
if version == 1,
    
    logLikmat_filt = zeros(nr_pat, 1);
    for pat = 1:nr_pat,
        
        logLikmat_filt(pat) = logLikmat(pat, I_pat(pat));;
        
    end
    logLik = sum(sum(logLikmat_filt));
    model_bestfit = atrophy_model(:, order_events(I_pat));
    
elseif version == 2,
    
    Likmat = zeros(nr_pat, nr_roi);
    for roi = 1:nr_roi,
        
        Likmat(:, roi) = prod(p_A_D(atrophy_model(:, order_events(roi)), :), 1).*...
            prod(p_noA_D(noatrophy_model(:, order_events(roi)), :), 1);
        
    end
    logLik = sum(log(sum(Likmat, 2)));
    if nargout > 1,
        
        logLikmat = zeros(nr_pat, nr_roi);
        for roi = 1:nr_roi,
            
            logLikmat(:, roi) = sum(log(p_A_D(atrophy_model(:, order_events(roi)), :)), 1) + ...
                sum(log(p_noA_D(noatrophy_model(:, order_events(roi)), :)), 1);
            
        end

        I_pat = zeros(nr_pat, 1);
        for pat = 1:nr_pat,
            
            I_pat_local = find(logLikmat(pat, :) == max(logLikmat(pat, :)));
            I_pat(pat) = I_pat_local(1);
            
        end
        model_bestfit = atrophy_model(:, order_events(I_pat));
        
    end
    
elseif version == 3,
    
    logLik_stage = zeros(nr_roi, 1);
    for stage = 1:nr_roi,
        
        
        logLik_stage(stage) = sum(prod(p_A_D(atrophy_model(:, stage), :), 1).*...
            prod(p_noA_D(noatrophy_model(:, stage), :), 1));
        
    end
    logLik = sum(log(logLik_stage));
    if nargout == 3,
        
        logLikmat = zeros(nr_pat, nr_roi);
        for roi = 1:nr_roi,
            
            for pat = 1:nr_pat,
                
                logLikmat(pat, roi) = logLikmat(pat, roi) + ...
                    sum(logp_A_D(atrophy_model(:, order_events(roi)), pat));
                logLikmat(pat, roi) = logLikmat(pat, roi) + ...
                    sum(logp_noA_D(noatrophy_model(:, order_events(roi)), pat));
                
            end
            
        end
        I_pat = zeros(nr_pat, 1);
        for pat = 1:nr_pat,
            
            I_pat_local = find(logLikmat(pat, :) == max(logLikmat(pat, :)));
            I_pat(pat) = I_pat_local(1);
            
        end
        model_bestfit = atrophy_model(:, order_events(I_pat));
        
    end
    
end      

