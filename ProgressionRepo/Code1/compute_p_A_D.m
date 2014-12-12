function [p_A_D, prob_pat, prob_con, span] = ...
    compute_p_A_D(atrophy_control, atrophy_patient, flag_mixt, flag_filt, ...
    flag_vis)

if ~exist('flag_vis', 'var'),
    
    flag_vis = 0;
    
end

span_min = min([min(atrophy_patient.data_mean, [], 2) ...
    min(atrophy_control.data_mean, [], 2)], [], 2);
span_max = max([max(atrophy_patient.data_mean, [], 2) ...
    max(atrophy_control.data_mean, [], 2)], [], 2);
span = [span_min span_max];

if flag_mixt == 1,

    atrophy_patient = atrophy_patient.data_mean;
    atrophy_control = atrophy_control.data_mean;
    [nr_roi, nr_pat] = size(atrophy_patient);
    p_A_D = zeros(nr_roi, nr_pat);

    Cmin = 1;
    Cmax = 5;
    it_vec = 100*ones(1, Cmax);
    for roi = 1:nr_roi,
        
        gmix_c = netlab_wrapper(atrophy_control(roi, :)', ...
            Cmin, Cmax, it_vec, 'full', 0);
        component_max = find(gmix_c.priors == max(gmix_c.priors));
        
        if flag_filt == 1,
            
            gmix_c_filt = gmm(1, 1, 'full');
            gmix_c_filt.centres = gmix_c.centres(component_max(1));
            gmix_c_filt.covars = gmix_c.covars(:, :, component_max(1));
            
        elseif flag_filt == 0,
            
            gmix_c_filt = gmix_c;
            
        end
        
        p_D_noA = gmmprob(gmix_c_filt, atrophy_patient(roi, :)');
        gmix_t = netlab_wrapper(cat(1, atrophy_control(roi, :)', ...
            atrophy_patient(roi, :)'), Cmin, Cmax, it_vec, 'full', 0);
        p_D = gmmprob(gmix_t, atrophy_patient(roi, :)');
        
        cdf_c = ksdensity(atrophy_control(roi, :), atrophy_patient(roi, :), ...
            'function', 'cdf');
        I_noatrophy = find(cdf_c > 0.5);
        p_D_noA(I_noatrophy) = max(cat(2, p_D_noA(I_noatrophy), ...
            p_D(I_noatrophy)), [], 2);
        p_D(I_noatrophy) = max(cat(2, p_D_noA(I_noatrophy), ...
            p_D(I_noatrophy)), [], 2);
        p_D_A = p_D - p_D_noA;
        p_A_D_roi = p_D_A./(p_D_A + p_D_noA);
        p_A_D_roi(p_A_D_roi < 0) = 0;
        p_A_D_roi(p_A_D_roi == 0) = eps;
        p_A_D(roi, :) = p_A_D_roi;
        fprintf('roi: %d in %d regions\n', roi, nr_roi);
        
    end
    
elseif flag_mixt == 0,
    
    atrophy_patient = atrophy_patient.data_mean;
    atrophy_control = atrophy_control.data_mean;
    [nr_roi, nr_pat] = size(atrophy_patient);
    p_A_D = zeros(nr_roi, nr_pat);
    for roi = 1:nr_roi,
        
        [mu_control(roi), std_control(roi)] = ...
            normfit(atrophy_control(roi, :));
        
    end
    for pat = 1:nr_pat,
        
        for roi = 1:nr_roi,
            
            [h, p] = ...
                ztest(atrophy_patient(roi, pat), ...
                mu_control(roi), std_control(roi), [], 'left');
            p_A_D(roi, pat) = p;
            
        end
        
    end
    thr1 = 0.05;
    I_thr = find(p_A_D > thr1);
    p_A_D = 1 - p_A_D;
    p_A_D(I_thr) = 0;
    prob_pat = p_A_D;
    prob_con = 1-p_A_D;
    
elseif flag_mixt == 2,
       
    [nr_roi, nr_pat] = size(atrophy_patient.data_mean);
    p_A_D = zeros(nr_roi, nr_pat);
    mu_control = zeros(nr_roi, 1);
    sig_control = zeros(nr_roi, 1);
    N_control = zeros(nr_roi, 1);
    for roi = 1:nr_roi,
        
        mu_local = atrophy_control.data_mean(roi, :);
        sig_local = atrophy_control.data_std(roi, :);
        N_local = atrophy_control.data_nr_vox(roi, :);
        mu_control(roi) = sum(N_local.*mu_local)/sum(N_local);
        sig_control(roi) = ...
            sqrt((sum(N_local.*((sig_local.^2) + (mu_local.^2)))/sum(N_local)) - ...
            (mu_control(roi).^2));
        N_control(roi) = sum(N_local);            
            
    end
    for pat = 1:nr_pat,
        
        for roi = 1:nr_roi,
            
            data_control = normrnd(mu_control(roi), sig_control(roi), N_control(roi), 1);
            data_patient = normrnd(atrophy_patient.data_mean(roi, pat), ...
                atrophy_patient.data_std(roi, pat), ...
                atrophy_patient.data_nr_vox(roi, pat), 1);
            [h, p] = ttest2(data_patient, data_control, [], 'left');                   
            p_A_D(roi, pat) = 1-p;
            
        end
        
    end

elseif flag_mixt == 3,
    
    atrophy_patient = atrophy_patient.data_mean;
    atrophy_control = atrophy_control.data_mean;
    [nr_roi, nr_pat] = size(atrophy_patient);
    nr_con = size(atrophy_control, 2);
    p_A_D = zeros(nr_roi, nr_pat);
    prob_pat = zeros(nr_roi, nr_pat);
    prob_con = zeros(nr_roi, nr_pat);

    Cmin = 1;
    Cmax = 5;
    it_vec = 100*ones(1, Cmax-Cmin);
    for roi = 1:nr_roi,
            
        data_control = atrophy_control(roi, :)';
        %         I_purge = find(data_control < 0.6);
        %         data_control(I_purge) = [];
        data_patient = atrophy_patient(roi, :)';
        %         I_purge = find(data_patient < 0.6);
        %         data_patient(I_purge) = [];
        data_tot = cat(1, data_control, data_patient);
        gmix_t = netlab_wrapper(data_tot, ...
            Cmin, Cmax, it_vec, 'full', 0);
        if gmix_t.ncentres == 1,
            
            [mu_control, std_control] = ...
                normfit(atrophy_control(roi, :));
            p_A_D_roi = zeros(nr_pat, 1);
            for pat = 1:nr_pat,
                
                [h, p] = ...
                    ztest(atrophy_patient(roi, pat), ...
                    mu_control, std_control, [], 'left');
                p_A_D_roi(pat) = p;
                
            end
            [thr1, thr2] = FDR(p_A_D_roi, 0.05);
            if isempty(thr1),
                
                thr1 = 0.001;
                
            end
            I_thr = find(p_A_D_roi > thr1);
            p_A_D_roi = 1-p_A_D_roi;
            p_A_D_roi(I_thr) = 0;
            p_A_D(roi, :) = p_A_D_roi;
            prob_pat(roi, :) = p_A_D_roi;
            prob_con(roi, :) = 1-p_A_D_roi;
            
        else

            post_c = gmmpost(gmix_t, atrophy_control(roi, :)');
            post_p = gmmpost(gmix_t, atrophy_patient(roi, :)');
            nr_controls_component = zeros(gmix_t.ncentres, 1);
            nr_patients_component = zeros(gmix_t.ncentres, 1);
            for c = 1:gmix_t.ncentres,
                
                I_clust = find(post_c(:, c) == max(post_c, [], 2));
                nr_controls_component(c) = length(I_clust)/nr_con;
                I_clust = find(post_p(:, c) == max(post_p, [], 2));
                nr_patients_component(c) = length(I_clust)/nr_pat;
                
            end
            
            component_control = find((nr_controls_component > 0.5*nr_patients_component));
            if isempty(component_control),
                
                component_control = find(nr_controls_component == max(nr_controls_component));
                
            end
            nr_component_control(roi) = length(component_control);
            component_atrophy = find(gmix_t.centres < ...
                min(gmix_t.centres(component_control)));
            activ_p = gmmactiv(gmix_t, data_patient);
            prob_comp = activ_p.*repmat(gmix_t.priors, size(activ_p, 1), 1);
            prob_comp(atrophy_patient(roi, :) >= ...
                min(gmix_t.centres(component_control)), component_atrophy) = 0;
            prob_pat_roi = sum(prob_comp(:, component_atrophy), 2);
            prob_pat(roi, :) = prob_pat_roi;
            prob_con(roi, :) = sum(prob_comp(:, component_control), 2);
            
            p_A_D_roi = sum(post_p(:, component_atrophy), 2);
            p_A_D_roi(atrophy_patient(roi, :) >= min(gmix_t.centres(component_control))) = 0;
            p_A_D(roi, :) = p_A_D_roi;
            
        end
        
        if flag_vis,
            
            [n_t, x_t] = hist(data_tot, 20);
            figure(1), clf
            subplot(1, 2, 1), hold on
            bar(x_t, n_t, 'k')
            p_c = gmmprob(gmix_t, atrophy_control(roi, :)');
            p_c = p_c*(max(n_t)/max(p_c));
            p_p = gmmprob(gmix_t, atrophy_patient(roi, :)');
            p_p = p_p*(max(n_t)/max(p_p));
            plot(atrophy_control(roi, :), p_c, '*g')
            plot(atrophy_patient(roi, :), p_p, '.r')
            axis([min(data_tot) max(data_tot) 0 max(n_t)])
            subplot(1, 2, 2), hold on
            plot(atrophy_patient(roi, :), p_A_D_roi, '.k')
            axis([min(data_tot) max(data_tot) 0 1])
            fprintf('roi: %d in %d regions\n', roi, nr_roi);
            pause
            
        else
            
            fprintf('roi: %d in %d regions\n', roi, nr_roi);
            
        end
        
    end
    nr_component_control
elseif flag_mixt == 4,
    
    atrophy_patient = atrophy_patient.data_mean;
    atrophy_control = atrophy_control.data_mean;
    [nr_roi, nr_pat] = size(atrophy_patient);
    p_A_D = zeros(nr_roi, nr_pat);

    Cmin = 1;
    Cmax = 5;
    it_vec = 10*ones(1, Cmax-Cmin);
    opts = foptions;
    opts(14) = 1e3;
    for roi = 1:nr_roi,
        
        data_control = atrophy_control(roi, :)';
        data_patient = atrophy_patient(roi, :)';
        data_tot = cat(1, data_control, data_patient);
        
        gmix_c = netlab_wrapper(data_control, 1, ...
            1, 10, 'full', 0);
        [gmix_t, BIC_t, gmix_struct] = ...
            gmmem_fixcomp_wrapper(gmix_c, data_tot, Cmax, max(it_vec), opts);
        post_c = gmmpost(gmix_t, data_control);
        post_p = gmmpost(gmix_t, data_patient);
        component_control = zeros(gmix_c.ncentres, 1);
        for c = 1:gmix_c.ncentres,
            
            component_control(c) = ...
                find(gmix_t.centres == gmix_c.centres(c));
            
        end
        component_atrophy = find(gmix_t.centres < ...
            min(gmix_t.centres(component_control)));
        p_A_D_roi = sum(post_p(:, component_atrophy), 2);
        p_A_D_roi(atrophy_patient(roi, :) >= min(gmix_t.centres(component_control))) = 0;
        p_A_D(roi, :) = p_A_D_roi;

        [n_c, x_c] = hist(data_control, 20);
        [n_p, x_p] = hist(data_patient, 20);
        figure(1), clf
        subplot(1, 2, 1), hold on
        bar(x_c, n_c, 'g')
        bar(x_p, n_p, 'r')
        p_c = gmmprob(gmix_t, atrophy_control(roi, :)');
        %         p_c = p_c*(max(n_t)/max(p_c));
        p_p = gmmprob(gmix_t, atrophy_patient(roi, :)');
        %         p_p = p_p*(max(n_t)/max(p_p));
        plot(atrophy_control(roi, :), p_c, '*g')
        plot(atrophy_patient(roi, :), p_p, '.r')
        %         axis([min(data_tot) max(data_tot) 0 max(n_t)])
        subplot(1, 2, 2), hold on
        plot(atrophy_patient(roi, :), p_A_D_roi, '.k')
        axis([min(data_tot) max(data_tot) 0 1])
        fprintf('roi: %d in %d regions\n', roi, nr_roi);
        pause
        
    end
    
elseif flag_mixt == 6,
    
    atrophy_patient = atrophy_patient.data_mean;
    atrophy_control = atrophy_control.data_mean;
    [nr_roi, nr_pat] = size(atrophy_patient);
    nr_con = size(atrophy_control, 2);
    p_A_D = zeros(nr_roi, nr_pat);
    prob_pat = zeros(nr_roi, nr_pat);
    prob_con = zeros(nr_roi, nr_pat);

    Cmin = 1;
    Cmax = 5;
    it_vec = 100*ones(1, Cmax-Cmin);
    for roi = 1:nr_roi,
            
        data_control = atrophy_control(roi, :)';
        %         I_purge = find(data_control < 0.6);
        %         data_control(I_purge) = [];
        data_patient = atrophy_patient(roi, :)';
        %         I_purge = find(data_patient < 0.6);
        %         data_patient(I_purge) = [];
        data_tot = cat(1, data_control, data_patient);
        gmix_t = netlab_wrapper(data_tot, ...
            Cmin, Cmax, it_vec, 'full', 0);
        if gmix_t.ncentres == 1,
            
            [mu_control, std_control] = ...
                normfit(atrophy_control(roi, :));
            p_A_D_roi = zeros(nr_pat, 1);
            for pat = 1:nr_pat,
                
                [h, p] = ...
                    ztest(atrophy_patient(roi, pat), ...
                    mu_control, std_control, [], 'left');
                p_A_D_roi(pat) = p;
                
            end
            [thr1, thr2] = FDR(p_A_D_roi, 0.05);
            if isempty(thr1),
                
                thr1 = 0.001;
                
            end
            I_thr = find(p_A_D_roi > thr1);
            p_A_D_roi = 1-p_A_D_roi;
            p_A_D_roi(I_thr) = 0;
            p_A_D(roi, :) = p_A_D_roi;
            prob_pat(roi, :) = p_A_D_roi;
            prob_con(roi, :) = 1-p_A_D_roi;
            
        else

            post_c = gmmpost(gmix_t, atrophy_control(roi, :)');
            post_p = gmmpost(gmix_t, atrophy_patient(roi, :)');
            nr_controls_component = zeros(gmix_t.ncentres, 1);
            nr_patients_component = zeros(gmix_t.ncentres, 1);
            for c = 1:gmix_t.ncentres,
                
                I_clust = find(post_c(:, c) == max(post_c, [], 2));
                nr_controls_component(c) = length(I_clust)/nr_con;
                I_clust = find(post_p(:, c) == max(post_p, [], 2));
                nr_patients_component(c) = length(I_clust)/nr_pat;
                
            end
            
            component_control = find(nr_controls_component == max(nr_controls_component));
            nr_component_control(roi) = length(component_control);
            component_atrophy = find(gmix_t.centres < ...
                min(gmix_t.centres(component_control)));
            activ_p = gmmactiv(gmix_t, data_patient);
            prob_comp = activ_p.*repmat(gmix_t.priors, size(activ_p, 1), 1);
            prob_comp(atrophy_patient(roi, :) >= ...
                min(gmix_t.centres(component_control)), component_atrophy) = 0;
            prob_pat_roi = sum(prob_comp(:, component_atrophy), 2);
            prob_pat(roi, :) = prob_pat_roi;
            prob_con(roi, :) = sum(prob_comp(:, component_control), 2);
            
            p_A_D_roi = sum(post_p(:, component_atrophy), 2);
            p_A_D_roi(atrophy_patient(roi, :) >= min(gmix_t.centres(component_control))) = 0;
            p_A_D(roi, :) = p_A_D_roi;
            
        end
        
        if flag_vis,
            
            [n_t, x_t] = hist(data_tot, 20);
            figure(1), clf
            subplot(1, 2, 1), hold on
            bar(x_t, n_t, 'k')
            p_c = gmmprob(gmix_t, atrophy_control(roi, :)');
            p_c = p_c*(max(n_t)/max(p_c));
            p_p = gmmprob(gmix_t, atrophy_patient(roi, :)');
            p_p = p_p*(max(n_t)/max(p_p));
            plot(atrophy_control(roi, :), p_c, '*g')
            plot(atrophy_patient(roi, :), p_p, '.r')
            axis([min(data_tot) max(data_tot) 0 max(n_t)])
            subplot(1, 2, 2), hold on
            plot(atrophy_patient(roi, :), p_A_D_roi, '.k')
            axis([min(data_tot) max(data_tot) 0 1])
            fprintf('roi: %d in %d regions\n', roi, nr_roi);
            pause
            
        else
            
            fprintf('roi: %d in %d regions\n', roi, nr_roi);
            
        end
        
    end
end

