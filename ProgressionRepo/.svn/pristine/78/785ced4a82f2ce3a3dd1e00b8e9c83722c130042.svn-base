function p_atrophy = compute_p_atrophy_mixt(atrophy_control, atrophy_patient)

[nr_roi, nr_control] = size(atrophy_control);
nr_patient = size(atrophy_patient, 2);
p_atrophy = zeros(nr_roi, nr_patient);
opts = foptions;
opts(14) = 1e4;
nr_it = 1e2;
for roi = 1:nr_roi,
    
    data_control = atrophy_control(roi, :);
    data_patient = atrophy_patient(roi, :);
    data_tot = [data_control data_patient];
    [n_c, x_c] = hist(data_control, 5);
    [n_p, x_p] = hist(data_patient, 5);
    [gmix_control{roi}, BIC_d, gmix_c_d] = ...
        netlab_wrapper(data_control', 1, 2, nr_it*ones(1, 2), 'full', 0);
    [gmix_tot_nonfix{roi}, BIC_nonfix, gmix_struct_nonfix] = ...
        netlab_wrapper(data_tot', 1, 2, nr_it*ones(1, 2), 'full', 0);
    [gmix_tot{roi}, BIC_extracomp, gmix_struct] = ...
        gmmem_fixcomp_wrapper(gmix_control{roi}, ...
        data_tot', 3, nr_it, opts);
    BIC_control2tot = gmmem_BIC(gmix_control{roi}, data_tot');
    if ((BIC_control2tot < min(BIC_extracomp)) || isempty(gmix_tot{roi})),
        
        p_atrophy(roi, :) = 0;
        
    else
        
        I_control = [1:gmix_control{roi}.ncentres]';
        I_nuis = find(gmix_tot{roi}.centres((gmix_control{roi}.ncentres+1):end) > ...
            max(gmix_tot{roi}.centres(I_control)));
        I_patient = 1:gmix_tot{roi}.ncentres;
        I_patient([I_control; I_nuis]) = [];
        if isempty(I_patient),
            
            p_atrophy(roi, :) = 0;
            
        else
            
            post = gmmpost(gmix_tot{roi}, data_patient');
            p_atrophy(roi, :) = sum(post(:, I_patient), 2);
            
        end
        
    end       
    fprintf('roi: %d in %d rois\n', roi, nr_roi)
    
end
