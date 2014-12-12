% [logLik, logLikmat, model_bestfit] = logLikelihoodAtrophyModel(p_A_D, ...
%     order_events, version)
% Version signifies the way likelihood is calculated
% version == 1: the stage with the highest likelihood for each subject is
% taken
% version == 2: The likelihood of all stages is summed for each subject,
% then the product over all subjects is taken
% version == 3 The likelihood for all patients is summed for each stage,
% then the product over all stages is taken

function [logLik, Likmat, model_bestfit] = ...
    logLikelihoodAtrophyModelCorrect(prob_pat, prob_con, ...
    order_events)

[nr_roi, nr_pat] = size(prob_pat);
event_mat = zeros(nr_roi, nr_roi);
for roi = 1:nr_roi,
    
    event_mat(order_events==roi, roi:end) = 1;
    
end
for k = 1:nr_roi,
    
    I_event{k} = find(event_mat(:, k) == 1);
    I_noevent{k} = find(event_mat(:, k) == 0);
    
end

Likmat = zeros(nr_roi, nr_pat);
for j = 1:nr_pat,
    
    for k = 1:nr_roi,
        
        Likmat(k, j) = prod(prob_pat(I_event{k}, j)).*prod(prob_con(I_noevent{k}, j));
        
    end
    
end
logLik = sum(log(sum(Likmat, 1)));

model_bestfit = zeros(nr_roi, nr_pat);
if nargout == 3,
    
    for j = 1:nr_pat,
        
        I = find(Likmat(:, j) == max(Likmat(:, j)));
        model_bestfit(:, j) = event_mat(:, I(1));
        
    end
    
end

