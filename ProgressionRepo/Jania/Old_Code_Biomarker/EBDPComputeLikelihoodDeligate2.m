% function to calculate likelihood of data given event or no event
% by fitting a mixture of Gaussians to the data.
% [likelihood, gmix_events] = ...
%     EBDPComputeLikelihood(data_patients, data_controls, version)
% version = 1: mixture model in which first a single Gaussian is fittied to
% the control distribution (in which no event has occurred)
% Then a mixture is fitted to the complete data, while keeping the no-event
% distribution fixed
%
% version = 2: free mixture model fitted to all of the data, some trickery
% to find the no-event and event components
% version = 4: mixture for Gaussian and uniform distribution
% version = 6: first fit a MoG to the control distributions, then fit a MoG
% to the patient distribution, with the control means and variances
% included and fixed. Each time new components are introduced with limited
% variance.
function [likelihood, gmix_struct] = ...
    EBDPComputeLikelihood(data_patients, data_controls, version)

[nr_pat, nr_roi] = size(data_patients);
nr_cont = size(data_controls, 1);
nr_clust_min = 1;
nr_clust_max = 5;
likelihood = zeros(nr_roi, size(data_patients, 1), 2);
opts = foptions;
opts(5) = 1;
opts(14) = 1e3;
opts(19) = 10; % nr_components_add_max
opts(20) = 0.6; % reduction_covars
opts(21) = 1e2; % #repetitions for adding components
data_tot = [data_controls; data_patients];
flag_direction_extracomponents = 1;
parm_mixtUnif.nr_it = 1e3;
I_patients = (nr_cont+1):(nr_cont+nr_pat);

if version == 5,
    
    p_sig = zeros(nr_roi, nr_pat);
    for roi = 1:nr_roi,
        
        [mu_control, std_control] = normfit(data_patients(:, roi));
        for pat = 1:nr_pat,
            
            [h, p] = ...
                ztest(data_patients(pat, roi), ...
                mu_control, std_control, [], 'left');
            p_sig(roi, pat) = p;
            
        end
        
    end
    %     [thr1, thr2] = FDR(p_sig(:), 0.05);
    thr1 = 0.05;
    likelihood_noevents = p_sig >= thr1;
    likelihood_events = p_sig < thr1;
    likelihood(:, :, 1) = likelihood_noevents;
    likelihood(:, :, 2) = likelihood_events;
    gmix_struct = [];
        
end          
            
for roi = 1:nr_roi,
    
    if version == 1,
        
        %     Fit Gaussian mixture
        gmix_controls = gmm(1, 1, 'full');
        gmix_controls = gmmem(gmix_controls, data_controls(:, roi), opts);
        gmix_tot = gmmem_fixcomp_wrapper(gmix_controls, data_tot(:, roi), ...
            2, 10, opts, flag_direction_extracomponents);
        gmix_struct.gmix_roi{roi} = gmix_tot;
        component_patients = find(gmix_tot.centres < gmix_controls.centres);
        if isempty(component_patients),
            
            [mu_control, std_control] = normfit(data_patients(:, roi));
            likelihood_local = normpdf(data_patients(:, roi), mu_control, std_control);
            p_sig = zeros(nr_pat, 1);
            for pat = 1:nr_pat,
                
                [h, p] = ...
                    ztest(data_patients(pat, roi), ...
                    mu_control, std_control, [], 'left');
                p_sig(pat) = p;
                
            end
            I_sig = find(p_sig < 0.01);
            I_nonsig = find(p_sig >= 0.01);
            likelihood(roi, I_nonsig, 1) = likelihood_local(I_nonsig);
            likelihood(roi, I_sig, 1) = 0;
            likelihood(roi, I_sig, 2) = 1./likelihood_local(I_sig);
            likelihood(roi, I_nonsig, 2) = 0;
            
        else
            
            gmix_patients = gmm(1, length(component_patients), 'full');
            gmix_patients.centres = gmix_tot.centres(component_patients);
            gmix_patients.covars = gmix_tot.covars(:, :, component_patients);
            gmix_patients.priors = gmix_tot.priors(component_patients);
            gmix_patients.priors = gmix_patients.priors/sum(gmix_patients.priors);
            likelihood(roi, :, 1) = gmmprob(gmix_controls, data_patients(:, roi));
            likelihood(roi, :, 2) = gmmprob(gmix_patients, data_patients(:, roi));
            
        end
        
    elseif version == 3,

        gmix_tot = netlab_wrapper(data_tot(:, roi), nr_clust_min, nr_clust_max, ...
            100*ones(1, nr_clust_max - nr_clust_min + 1), 'full', 0);
        gmix_struct.gmix_roi{roi} = gmix_tot;
        post_controls = gmmpost(gmix_tot, data_controls(:, roi));
        nr_controls_component = zeros(gmix_tot.ncentres, 1);
        for c = 1:gmix_tot.ncentres,
            
            I_c = find(post_controls(:, c) == max(post_controls, [], 2));
            nr_controls_component(c) = length(I_c);
            
        end
        component_controls = find(nr_controls_component == ...
            max(nr_controls_component));
        component_controls = component_controls(1);
        component_patients = find(gmix_tot.centres < ...
            gmix_tot.centres(component_controls));
        if isempty(component_patients), % if we don't really find a patient component, 
            % I'm reverting to a z-test on the control disbtribution
                        
            [mu_control, std_control] = ...
                normfit(data_controls(:, roi));
            likelihood_local = normpdf(data_patients(:, roi), mu_control, std_control);
            p_sig = zeros(nr_pat, 1);
            for pat = 1:nr_pat,
                
                [h, p] = ...
                    ztest(data_patients(pat, roi), ...
                    mu_control, std_control, [], 'left');
                p_sig(pat) = p;
                
            end
            [thr1, thr2] = FDR(p_sig, 0.05); % False discovery rate correction 
            % for multiple comparisons, which is just a way to avoid arbitrary 
            % thresholds (other then deciding that p = 0.05 FDR-corrected is 
            % a good threshold
            if isempty(thr1),
                
                thr1 = 0.001; % if this fails, set an arbitrary threshold
                
            end
            
            
            % There's no really good way to set the probability of the
            % measurement, given that the event has occurred, so I'm just
            % using some arbitrary values here. And yes, these values
            % influence where the event will end up... 
            
            I_sig = find(p_sig < thr1);
            I_nonsig = find(p_sig >= thr1);
            likelihood(roi, I_nonsig, 1) = likelihood_local(I_nonsig);
            likelihood(roi, I_sig, 1) = 0.001;
            likelihood(roi, I_sig, 2) = 1;
            likelihood(roi, I_nonsig, 2) = 0.001;
            
        else
            
            gmix_controls = gmm(1, 1, 'full');
            gmix_controls.centres = gmix_tot.centres(component_controls);
            gmix_controls.covars = gmix_tot.covars(:, :, component_controls);
            
            gmix_patients = gmm(1, length(component_patients), 'full');
            gmix_patients.centres = gmix_tot.centres(component_patients);
            gmix_patients.covars = gmix_tot.covars(:, :, component_patients);
            gmix_patients.priors = gmix_tot.priors(component_patients);
            gmix_patients.priors = gmix_patients.priors/sum(gmix_patients.priors);
            
            likelihood(roi, :, 1) = gmmprob(gmix_controls, data_patients(:, roi));
            likelihood(roi, :, 2) = gmmprob(gmix_patients, data_patients(:, roi));           
            
        end
        
    elseif version == 2,

        gmix_tot = netlab_wrapper(data_tot(:, roi), nr_clust_min, nr_clust_max, ...
            100*ones(1, 5), 'full', 0);
        gmix_struct.gmix_roi{roi} = gmix_tot;
        post_controls = gmmpost(gmix_tot, data_controls(:, roi));
        nr_controls_component = zeros(gmix_tot.ncentres, 1);
        for c = 1:gmix_tot.ncentres,
            
            I_c = find(post_controls(:, c) == max(post_controls, [], 2));
            nr_controls_component(c) = length(I_c);
            
        end
        component_controls = find(nr_controls_component == ...
            max(nr_controls_component));
        component_controls = component_controls(1);
        component_patients = find(gmix_tot.centres < ...
            gmix_tot.centres(component_controls));
        if isempty(component_patients), % if we don't really find a patient component, 
            % I'm reverting to a z-test on the control disbtribution
                        
            [mu_control, std_control] = ...
                normfit(data_controls(:, roi));
            likelihood_local = normpdf(data_patients(:, roi), mu_control, std_control);
            p_sig = zeros(nr_pat, 1);
            for pat = 1:nr_pat,
                
                [h, p] = ...
                    ztest(data_patients(pat, roi), ...
                    mu_control, std_control, [], 'left');
                p_sig(pat) = p;
                
            end
            % There's no really good way to set the probability of the
            % measurement, given that the event has occurred, so I'm just
            % using some arbitrary values here. And yes, these values
            % influence where the event will end up... 
            
            I_sig = find(p_sig < 0.01);
            I_nonsig = find(p_sig >= 0.01);
            likelihood(roi, I_nonsig, 1) = likelihood_local(I_nonsig);
            likelihood(roi, I_sig, 1) = 0.1;
            likelihood(roi, I_sig, 2) = 5;
            likelihood(roi, I_nonsig, 2) = 0.1;
            
        else
            
            gmix_controls = gmm(1, 1, 'full');
            gmix_controls.centres = gmix_tot.centres(component_controls);
            gmix_controls.covars = gmix_tot.covars(:, :, component_controls);
            
            gmix_patients = gmm(1, length(component_patients), 'full');
            gmix_patients.centres = gmix_tot.centres(component_patients);
            gmix_patients.covars = gmix_tot.covars(:, :, component_patients);
            gmix_patients.priors = gmix_tot.priors(component_patients);
            gmix_patients.priors = gmix_patients.priors/sum(gmix_patients.priors);
            
            likelihood(roi, :, 1) = gmmprob(gmix_controls, data_patients(:, roi));
            likelihood(roi, :, 2) = gmmprob(gmix_patients, data_patients(:, roi));           
            I_con = find(data_patients(:, roi) > gmix_controls.centres);
            likelihood(roi, I_con, 1) = max(squeeze(likelihood(roi, :, 1)));
            likelihood(roi, I_con, 2) = min(squeeze(likelihood(roi, :, 2)));
            
        end
        
    elseif version == 4,

        gmix_controls = gmm(1, 1, 'full');
        gmix_controls = gmmem(gmix_controls, data_controls(:, roi), opts);
        parm_mixtUnif.mu_init = gmix_controls.centres;
        parm_mixtUnif.sig2_init = gmix_controls.covars;
        % Here I'm setting the range of the uniform distribution. I want it
        % to run roughly from the mean of the control distribution (making
        % the likelihood that an event has happened given atrophy values
        % lower than this mean zero) and some extreme value. I'm setting it
        % to 0.7, which would mean a 30 % decrease in volume. This value
        % unfortunately has an effect on the likelihood estimation, so any
        % thoughts about better ways of setting this variable would be
        % appreciated...
        parm_mixtUnif.a = min(data_tot(:, roi)); 
        parm_mixtUnif.b = gmix_controls.centres;
        parm = fitFixedGaussUnif(data_tot(:, roi), parm_mixtUnif);
        parm.data = data_tot(:, roi);
        gmix_struct.gmix_roi{roi} = parm;
        
        likelihood(roi, :, 1) = parm.prob_N(I_patients);
        likelihood(roi, :, 2) = parm.prob_U(I_patients);
            
       
    end
    if version == 6,
        
        [gmix_controls, BIC_controls] = netlab_wrapper(data_controls(:, roi), ...
            1, 5, 100*ones(5, 1), 'spherical', 0);
        gmix_patients = gmmem_patient_components(gmix_controls, ...
            data_patients(:, roi), opts);
        
        gmix_noevent = gmix_controls;
        nr_components_add = gmix_patients.ncentres - gmix_controls.ncentres;
        gmix_event = gmm(1, nr_components_add, 'spherical');
        gmix_event.centres = ...
            gmix_patients.centres((gmix_controls.ncentres+1):end);
        gmix_event.covars = ...
            gmix_patients.covars((gmix_controls.ncentres+1):end);
        gmix_event.priors = ...
            gmix_patients.priors((gmix_controls.ncentres+1):end);
        gmix_event.priors = gmix_event.priors/sum(gmix_event.priors);
        likelihood(roi, :, 1) = gmmprob(gmix_noevent, data_patients(:, roi));
        likelihood(roi, :, 2) = gmmprob(gmix_event, data_patients(:, roi));
        
        gmix_struct.gmix_controls{roi} = gmix_controls;
        gmix_struct.gmix_patients{roi} = gmix_patients;
        
    end     
    fprintf('roi: %d in %d rois\n', roi, nr_roi)
    
    if version == 7
        gmix_controls = gmm(1, 1, 'full');
        gmix_controls = gmmem(gmix_controls, data_controls(:, roi), opts);
        
        gmix_patients = gmm(1, 1, 'full');
        gmix_patients = gmmem(gmix_patients, data_patients(:, roi), opts);
        
        likelihood(roi, :, 1) = gmmprob(gmix_controls, data_patients(:, roi));
        likelihood(roi, :, 2) = gmmprob(gmix_patients, data_patients(:, roi));
        
        gmix_struct.gmix_controls{roi} = gmix_controls;
        gmix_struct.gmix_patients{roi} = gmix_patients;
        
        
        
    end
    
    if version==8
         gmix_controls = gmm(1, 1, 'full');
         gmix_controls = gmmem(gmix_controls, data_controls(:, roi), opts);

         
         [gmix_tot, ~, gmix_struct_dummy]= gmmem_fixcomp_wrapper(gmix_controls, data_tot(:, roi), ...
            1, 10, opts, flag_direction_extracomponents);
        
        if gmix_tot.ncentres<2
            gmix_tot=gmix_struct_dummy{2};
            fprintf('%s','original BIC gmix had only one component');
        end
        
        gmix_patients = gmm(1, 1, 'full');
        gmix_patients.centres = gmix_tot.centres(2);
        gmix_patients.covars = gmix_tot.covars(:, :, 2);
        gmix_patients.priors = 1;

             
        likelihood(roi, :, 1) = gmmprob(gmix_controls, data_patients(:, roi));
        %likelihood(roi, :, 2) = gmmprob(gmix_tot, data_patients(:, roi));
        likelihood(roi, :, 2) = gmmprob(gmix_patients, data_patients(:, roi));
        
        gmix_struct.gmix_controls{roi} = gmix_controls;
        gmix_struct.gmix_patients{roi} = gmix_patients;
        gmix_struct.gmix_tot{roi}=gmix_tot;

        %% set a minimum likelihood
%         like=likelihood(roi,:,:);
%         like(like<1e-3)=1e-3;
%         likelihood(roi,:,:)=like;

        %set a minimum accepted likelihood equivalent to fitting a uniform
        %distribution with a fixed likelihood for all data
%          like1=likelihood(roi, :, 1);
%          like2=likelihood(roi, :, 2);
%          indx1=find(like1<1e-3);
%          indx2=find(like2<1e-3);
%          goodIndx1=setdiff([1:nr_pat],indx1);
%          goodIndx2=setdiff([1:nr_pat],indx2);
%          min_acc_like_patient1=min(like1(goodIndx1));
%          min_acc_like_patient2=min(like2(goodIndx2));
%          
%         
%          like1(like1<min_acc_like_patient1)=min_acc_like_patient1;
%          likelihood(roi, :, 1)=like1;
%          like2(like2<min_acc_like_patient2)=min_acc_like_patient2;
%          likelihood(roi, :, 2)=like2; 
%          
       %% force the datapoints in rih=ght tail to have higher Pr(not event) or equal Pr(not event) and Pr(event)
        tail_indx=find(data_patients(:,roi)>gmix_controls.centres);
        ratio=likelihood(roi,tail_indx,2)./likelihood(roi,tail_indx,1);
        indx2=find(ratio>1);
        likelihood(roi,tail_indx(indx2),1)=likelihood(roi,tail_indx(indx2),2);
         
%         %% set a maximumn likelilhood ratio
%         T=[2.2];
%         ratio=abs(log(likelihood(roi, :, 1))-log(likelihood(roi, :, 2)));
%         indx=find(ratio>log(T));
%          
%         for cIndx=1:length(indx)
%             thisLike=likelihood(roi,indx(cIndx),:);
%             %if you are in the far left of the distributrion
%             if data_patients(indx(cIndx))<gmix_controls.centres
%                 minIndx=find(thisLike==min(thisLike));
%                 mxIndx=find(thisLike==max(thisLike));
%                 thisLike(minIndx)=thisLike(mxIndx)/T;
%             else
%             %if you are in the far right of the distribution
%                 minIndx=find(thisLike==min(thisLike));
%                 mxIndx=find(thisLike==max(thisLike));
%                 thisLike(mxIndx)=thisLike(minIndx)*T; 
%             end
%               likelihood(roi,indx(cIndx),:)=thisLike;
%         end
%         
        
        
%        %% assign a maximum accepted likelihood
%                  
%          max_acc_like_p=gmmprob(gmix_patients,gmix_patients.centres);
%          max_acc_like_c=gmmprob(gmix_controls,gmix_controls.centres);
%         
%          max_acc_like=[min([max_acc_like_c,max_acc_like_p])];
%          
%          like=likelihood(roi, :, :);
%          like(like>max_acc_like)=max_acc_like;
%          likelihood(roi, :, :)=like;
        
%         like=likelihood(roi,:,:);
%         like(like<1e-3)=1e-3;
%         likelihood(roi,:,:)=like;

    end
    
    if version==9
         gmix_controls = gmm(1, 1, 'full');
         gmix_controls = gmmem(gmix_controls, data_controls(:, roi), opts);
         
         [gmix_tot, ~, gmix_struct_dummy]= gmmem_fixcomp_wrapper(gmix_controls, data_tot(:, roi), ...
            1, 10, opts, flag_direction_extracomponents);
        
        if gmix_tot.ncentres<2
            gmix_tot=gmix_struct_dummy{2};
            fprintf('%s','original BIC gmix had only one component');
        end
        
        gmix_patients = gmm(1, 1, 'full');
        gmix_patients.centres = gmix_tot.centres(2);
        gmix_patients.covars = gmix_tot.covars(:, :, 2);
        gmix_patients.priors = 1;
        
        likelihood(roi, :, 1) = gmmprob(gmix_controls, data_patients(:, roi));
        likelihood(roi, :, 2) = gmmprob(gmix_tot, data_patients(:, roi));
        
        gmix_struct.gmix_controls{roi} = gmix_controls;
        gmix_struct.gmix_patients{roi} = gmix_patients;
        gmix_struct.gmix_tot{roi}=gmix_tot;
        
    end
    
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function gmix_patients = gmmem_patient_components(gmix_controls, ...
    data_patients, opts)

nr_components_add_max = opts(19);
covars_controls = gmix_controls.covars;
covars_patients_max = opts(20)*mean(covars_controls);
BIC_components_add = zeros(nr_components_add_max, 1);
for nr_components_add = 1:nr_components_add_max,
    
    gmix_add{nr_components_add} = gmmem_add_components(gmix_controls, data_patients, ...
        nr_components_add, opts, covars_patients_max);
    BIC_components_add(nr_components_add) = ...
        gmmem_BIC(gmix_add{nr_components_add}, data_patients);
    
end
nr_components_opt = find(BIC_components_add == min(BIC_components_add));
gmix_patients = gmix_add{nr_components_opt};
    

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function gmix_add = gmmem_add_components(mix_fix, x, nr_components_add, ...
    opts, covars_patients)

[ndata, ndim] = size(x);

nr_rep = opts(21);
range_x = max(x) - min(x);
min_x = min(x);

BIC_rep = zeros(nr_rep, 1);
for it_rep = 1:nr_rep,
    
    mix = gmm(1, mix_fix.ncentres + nr_components_add, 'spherical');
    mix.centres(1:mix_fix.ncentres) = mix_fix.centres;
    mix.centres((mix_fix.ncentres + 1):end) = ...
        range_x*rand(nr_components_add, 1) + min_x;
    mix.covars(1:mix_fix.ncentres) = mix_fix.covars;
    mix.covars((mix_fix.ncentres + 1):end) = 0.2*covars_patients;
    mix.priors(1:mix_fix.ncentres) = mix_fix.priors;
    mix.priors = mix.priors/sum(mix.priors);

    % Sort out the options
    niters = opts(14);
    
    % Main loop of algorithm
    n = 1;
    flg_continue = 1;
    while (n <= niters) && flg_continue,
        
        % Calculate posteriors based on old parameters
        [post, act] = gmmpost(mix, x);
        
        prob = act*(mix.priors)';
        e = - sum(log(prob));
        if (n > 1 && abs(e - eold) < opts(3))
            opts(8) = e;
            flg_continue = 0;
        else
            eold = e;
        end
        
        % Adjust the new estimates for the parameters
        new_pr = sum(post, 1);
        new_c = post' * x;
  
        % Now move new estimates to old parameter vectors
        mix.priors = new_pr ./ ndata;
  
        mix.centres = new_c ./ (new_pr' * ones(1, mix.nin));
        mix.centres(1:mix_fix.ncentres) = mix_fix.centres;
        n2 = dist2(x, mix.centres);
        for j = 1:mix.ncentres
            v(j) = (post(:,j)'*n2(:,j));
        end
        mix.covars = ((v./new_pr))./mix.nin;
        mix.covars(1:mix_fix.ncentres) = mix_fix.covars;
        for j = (mix_fix.ncentres+1):mix.ncentres,
            
            if mix.covars(j) > covars_patients,
                
                mix.covars(j) = covars_patients;
                
            end
            
            if mix.covars(j) < 1e-5,
                
                mix.covars(j) = 1;
                
            end
            
        end  
        if find(isnan(mix.centres) + isnan(mix.covars')),
            
            mix = gmm(1, mix_fix.ncentres + nr_components_add, 'spherical');
            mix.centres(1:mix_fix.ncentres) = mix_fix.centres;
            mix.centres((mix_fix.ncentres + 1):end) = ...
                range_x*rand(nr_components_add, 1) + min_x;
            mix.covars(1:mix_fix.ncentres) = mix_fix.covars;
            mix.covars((mix_fix.ncentres + 1):end) = 0.1*covars_patients;
            mix.priors(1:mix_fix.ncentres) = mix_fix.priors;
            mix.priors = mix.priors/sum(mix.priors);
            
        end
            
        n = n + 1;
    end
    gmix_rep{it_rep} = mix;
    BIC_rep(it_rep) = gmmem_BIC(gmix_rep{it_rep}, x);
    
end
it_opt = find(BIC_rep == min(BIC_rep));
gmix_add = gmix_rep{it_opt};