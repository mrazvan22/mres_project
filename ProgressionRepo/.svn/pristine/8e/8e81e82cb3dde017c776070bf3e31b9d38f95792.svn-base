function [est_event_order posteriorCase marginalPost] = EBDPMCMCTestWithParams1kCorrectClosedForm(data_controls, data_patients)

[nr_pat, nr_events] = size(data_patients);

flag_sum = 2;
flag_visualize = 0;
nr_it_params=1e4;

S=perms([1:nr_events]);
nr_cases=size(S,1);

log_probXGivenCase=zeros(nr_cases,nr_it_params);
thisMean_noevent_new =zeros(nr_events,nr_it_params);
thisMean_event_new=zeros(nr_events,nr_it_params);
thisCov_noevent_new=zeros(nr_events,nr_it_params);
thisCov_event_new=zeros(nr_events,nr_it_params);

prior_thisMean_noevent_new=zeros(nr_events,nr_it_params);
prior_thisMean_event_new=zeros(nr_events,nr_it_params);
prior_thisCov_noevent_new=zeros(nr_events,nr_it_params);
prior_thisCov_event_new=zeros(nr_events,nr_it_params);

likelihood_events_new=zeros(nr_events,nr_pat,2,nr_it_params);


%% First find the maximum likelihood parameters

% Get an initial estimate of the means and variances from the data
version_likelihood = 13;  %version_likelihood==13 fit a Mixture of Gaussians several
                          %times with different random initializations
nAttempts=5;              %number of times to fit MoG
threshold_flag=0;
[likelihood_events_est, gmix_struct_est] = ...
    EBDPComputeLikelihood(data_patients, data_controls, version_likelihood,threshold_flag,nAttempts);


%%  Evaluate the priors for current Mu and sigma
for roi=1:nr_events
    
    %% set means using mixGaussFit
    thisMean_noevent(roi)=gmix_struct_est.gmix_controls{roi}.mean;
    thisMean_event(roi)=gmix_struct_est.gmix_patients{roi}.mean;
    thisCov_noevent(roi)=gmix_struct_est.gmix_controls{roi}.cov;
    thisCov_event(roi)=gmix_struct_est.gmix_patients{roi}.cov;
    
    %% set the hyperparameters of the priors over Mu and Sigma
    %set the prior normal distribution centered over current mean
    init_mean_noevent(roi)=thisMean_noevent(roi);
    init_mean_event(roi)=thisMean_event(roi);
    init_cov_noevent(roi)=thisCov_noevent(roi);
    init_cov_event(roi)=thisCov_event(roi);
    
    
    %% set the sigma for event and noevent based on minimum and maximum
    %  value of the data_atient and data_control distributions separately
    
    %set the sigma of the prior normal distribution such that
    %min and max values of current biomarker are within + / - 3sigma
    minC=min(data_controls(:,roi));
    maxC=max(data_controls(:,roi));
    sigmaMin=(init_mean_noevent(roi)-minC)/3;
    sigmaMax=(maxC-init_mean_noevent(roi))/3;
    prior_sigma_noevent(roi)=max([sigmaMin,sigmaMax]);
    
    minP=min(data_patients(:,roi));
    maxP=max(data_patients(:,roi));
    sigmaMin=(init_mean_event(roi)-minP)/3;
    sigmaMax=(maxP-init_mean_event(roi))/3;
    prior_sigma_event(roi)=max([sigmaMin,sigmaMax]);
    
    
    %% Set the prior over mean to be Gaussian and the prior over variance to be a unifrom prior
    prior_thisMean_noevent(roi)=normpdf(thisMean_noevent(roi),init_mean_noevent(roi),prior_sigma_noevent(roi));
    prior_thisMean_event(roi)=normpdf(thisMean_event(roi),init_mean_event(roi),prior_sigma_event(roi));
    prior_thisCov_noevent(roi)=unifpdf(thisCov_noevent(roi),0,1e10);
    prior_thisCov_event(roi)=unifpdf(thisCov_event(roi),0,1e10);
    
    %% Visualize the likelihood and the prior distributions
    if flag_visualize
         
        %close all
        figure;
        %plot for no-event
        like_noevent_dist=normpdf(data_controls(:,roi),thisMean_noevent(roi),sqrt(thisCov_noevent(roi)));
        xmu0=linspace(min(data_controls(:,roi)),max(data_controls(:,roi)));
        prior_dist_mean_noevent=normpdf(xmu0,init_mean_noevent(roi),prior_sigma_noevent(roi));
        xsig0=[0.1:0.1:10];
        prior_dist_sigma_noevent=unifpdf(xsig0,0,realmax);
        cum_dist_prior_sigma_noevent=unifcdf(xsig0,0,realmax);
    
    
        figure
        subplot(1,4,1);plot(data_controls(:,roi),like_noevent_dist,'g.');title('Likelihood Dist No-event');
        hold on
        subplot(1,4,2);plot(xmu0,prior_dist_mean_noevent,'b.');title('Prior Dist over Mean')
        subplot(1,4,3);plot(xsig0,prior_dist_sigma_noevent,'m.');title('Prior Dist over Sigma');
        subplot(1,4,4);plot(xsig0,cum_dist_prior_sigma_noevent,'m.-');title('Cumulative Dist Prior Sigma');
        set(gcf,'Color',[ 1 1 1]);
    
        %plot for event
        figure
        like_event_dist=normpdf(data_patients(:,roi),thisMean_event(roi),sqrt(thisCov_event(roi)));
        xmu0=linspace(min(data_patients(:,roi)),max(data_patients(:,roi)));
        prior_dist_mean_event=normpdf(xmu0,init_mean_event(roi),prior_sigma_event(roi));
        xsig0=[0.1:0.1:10];
        prior_dist_sigma_event=unifpdf(xsig0,0,realmax);
        cum_dist_prior_sigma_event=unifcdf(xsig0,0,realmax);
    
        subplot(1,4,1);plot(data_patients(:,roi),like_event_dist,'r.');title('Likelihood Dist Event');
        hold on
        subplot(1,4,2);plot(xmu0,prior_dist_mean_event,'b.');title('Prior Dist over Mean')
        subplot(1,4,3);plot(xsig0,prior_dist_sigma_event,'m.');title('Prior Dist over Sigma');
        subplot(1,4,4);plot(xsig0,cum_dist_prior_sigma_event,'m.-');title('Cumulative Dist Prior Sigma');
        set(gcf,'Color',[ 1 1 1]);
    
    
        % compar the estimated dis vs true dist
        figure
        p_controls = like_noevent_dist;
        p_patients = like_event_dist;
        [hist_c, x_c] = ksdensity(data_controls(:, roi));
        [hist_p, x_p] = ksdensity(data_patients(:, roi));
    
        subplot(121), hold on
        plot(x_c, hist_c, 'g');
        plot(x_p, hist_p, 'r');
    
        subplot(122), hold on
        plot(data_controls(:, roi), p_controls, 'g.')
        plot(data_patients(:, roi), p_patients, 'r.');
        set(gcf,'Color',[ 1 1 1]);
    end
end

%% Find the likelihood of data given each ordering for the max like params

%Note I replace the log_likelihood_current of maximum likelihood with that
%of when also including the prior over parameters
for roi=1:nr_events
    %find the current likelihood event no event
    likelihood_events(roi, :, 1) = getGaussProb(data_patients(:, roi)',thisMean_noevent(roi),sqrt(thisCov_noevent(roi)));
    likelihood_events(roi, :, 2) = getGaussProb(data_patients(:, roi)',thisMean_event(roi),sqrt(thisCov_event(roi)));
end

% for each ordering
for c_case=1:nr_cases
    this_event_order=S(c_case,:);
    log_probXGivenCase(c_case,1) = findLikelihoodCaseGivenParams(likelihood_events, ...
        this_event_order, flag_sum, prior_thisMean_noevent, ...
        prior_thisMean_event,prior_thisCov_noevent,prior_thisCov_event);
    
end


%% Generate a list of possible Mu s and sigmas and the proposal distribution

%set the step size
change_step_muNE=0.02;
change_step_muE=0.02;
change_step_sigNE=0.02;
change_step_sigE=0.02;


min_sigma=0.09;

%  calculate a proposal distribution in this case based on the variance of the total data
for roi=1:nr_events
    sigma_event(roi)=std(data_patients(:,roi))^2;
    if sigma_event(roi)<min_sigma
        sigma_event(roi)=min_sigma;
    end
    sigma_noevent(roi)=std(data_controls(:,roi))^2;
    if sigma_noevent(roi)<min_sigma
        sigma_noevent(roi)=min_sigma;
    end
end


%initialize current parameters
for roi=1:nr_events
    thisMean_noevent_current(roi)=thisMean_noevent(roi);
    thisMean_event_current(roi)=thisMean_event(roi);
    thisCov_noevent_current(roi)= thisCov_noevent(roi);
    thisCov_event_current(roi)= thisCov_event(roi);
end

%% find the likelihood of data given each order and the new sampled parameters
for cParams = 2:nr_it_params
    
    if ~mod(cParams,100)
        cParams
    end
    % varying the means and the variances
    for roi=1:nr_events       
        
        %Update Samples from already accepted parameters    
        thisMean_noevent_new(roi,cParams)=thisMean_noevent_current(roi)+ sigma_noevent(roi)*change_step_muNE*randn;
        thisMean_event_new(roi,cParams)=thisMean_event_current(roi)+ sigma_event(roi)*change_step_muE*randn;
        thisCov_noevent_new(roi,cParams)=thisCov_noevent_current(roi)+ sigma_noevent(roi)*change_step_sigNE*randn;
        thisCov_event_new(roi,cParams)=thisCov_event_current(roi)+ sigma_event(roi)*change_step_sigE*randn;
        
        prior_thisMean_noevent_new(roi,cParams)=normpdf(thisMean_noevent_new(roi,cParams),init_mean_noevent(roi),prior_sigma_noevent(roi));
        prior_thisMean_event_new(roi,cParams)=normpdf(thisMean_event_new(roi,cParams),init_mean_event(roi),prior_sigma_event(roi));
        prior_thisCov_noevent_new(roi,cParams)=unifpdf(thisCov_noevent_new(roi,cParams),0,1e10);
        prior_thisCov_event_new(roi,cParams)=unifpdf(thisCov_event_new(roi,cParams),0,1e10);
               
        % constrain the Mean Event to always be smaller than Mean No-event
        if thisMean_event_new(roi,cParams)>thisMean_noevent_new(roi,cParams)
            prior_thisMean_noevent_new(roi,cParams)=0;
            prior_thisMean_event_new(roi,cParams)=0;
        end
             
        %find the current likelihood event no event
        likelihood_events_new(roi, :, 1,cParams) = getGaussProb(data_patients(:, roi)',thisMean_noevent_new(roi,cParams),sqrt(thisCov_noevent_new(roi,cParams)));
        likelihood_events_new(roi, :, 2,cParams) = getGaussProb(data_patients(:, roi)',thisMean_event_new(roi,cParams),sqrt(thisCov_event_new(roi,cParams)));
        
    end
    
    thisLikeEventsNew=squeeze(likelihood_events_new(:,:,:,cParams));
    thisPriorMeanNoEventNew=prior_thisMean_noevent_new(:,cParams);
    thisPriorMeanEventNew=prior_thisMean_event_new(:,cParams);
    thisPriorCovNoEventNew=prior_thisCov_noevent_new(:,cParams);
    thisPriorCovEventNew=prior_thisCov_event_new(:,cParams);
    
    % for each ordering
    for c_case=1:nr_cases
        
        this_event_order_new=S(c_case,:);
        log_probXGivenCase(c_case,cParams) = findLikelihoodCaseGivenParams(thisLikeEventsNew, ...
            this_event_order_new, flag_sum, thisPriorMeanNoEventNew', ...
            thisPriorMeanEventNew',thisPriorCovNoEventNew',thisPriorCovEventNew');
    end
    
end % end cParams iter

%calculate the posterior probability given each ordering and parameter

mx=max(log_probXGivenCase(:));
log_probXGivenCase_scaled=log_probXGivenCase-mx;
posteriorCase=exp(log_probXGivenCase_scaled)/sum(exp(log_probXGivenCase_scaled(:)));

%% Visualize the results
figure;
[X,Y] = meshgrid(1:100:1e4,1:nr_cases);
%draw the posterior probability of ordering and params
surf(X,Y,posteriorCase(:,1:100:end));
xlabel('Index of the Parameters Mu and Sigma');
ylabel('Ordering Index');
title('Posterior probability of Ordering Given Each Parameter');
set(gcf, 'Color',[1 1 1]);
%select the ordering with the highest posterior probability
[i,j]=find(posteriorCase==max(posteriorCase(:)));
est_event_order=S(i,:);

%find the marginal posterior distribution (i.e. integrate over the parameters Mu and Sigma)
marginalPost=sum(posteriorCase,2);

%plot the top 10% most likeliy orderings
top_n=max(1,round(nr_cases*0.1));
[top_n_order_indx]=findNmostLikelyOrderings(marginalPost,top_n);
most_likely_ordering=S(top_n_order_indx,:);
for c_event=1:nr_events
    for c_position=1:nr_events
        pos_var(c_event,c_position)=sum(most_likely_ordering(:,c_position)==est_event_order(c_event));
    end
end
figure
%plot the posterior probability of ordering marginalized over the parameters
plot(marginalPost);
set(gcf, 'Color',[1 1 1]);
xlabel('Ordering Index');
ylabel('Posterior Probability');
title('Posterior probability over ordering marginalized over parameters');

%plot the positional variance of the ordering
figure
imagesc(pos_var);
map_inverted = repmat(linspace(1, 0, 64)', 1, 3);
colormap(map_inverted)
axis square
set(gca, ...
    'YTick', [1:nr_events], 'YTickLabel',est_event_order, ...
    'YGrid', 'on')
hXLabel = xlabel('model stage');
hYLabel = ylabel('region');
set([hXLabel hYLabel], ...
    'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'demi');
set(gcf, 'Color',[1 1 1]);
title('Positional Variance')
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [log_likelihood_tot, class_pat] = findLikelihoodCaseGivenParams(likelihood_events, ...
    event_order, flag_sum, prior_thisMean_noevent, ...
    prior_thisMean_event,prior_thisCov_noevent,prior_thisCov_event)

[nr_events, nr_pat] = size(likelihood_events(:, :, 1));

likelihood_nonclinevents = likelihood_events;
likelihood_nonclinevents(likelihood_nonclinevents == 0) = ...
    min(likelihood_nonclinevents(likelihood_nonclinevents ~= 0));
log_likelihood_nonclinevents = log(likelihood_nonclinevents);

nr_phase = nr_events;
event_pattern = zeros(nr_events, nr_events);

%% *** NOTE **** here event pattern represents the given ordering
% e.g. if event_order= [1 3 2] this means event 1 happens first, event 3
% happend second and event 2 happend last
for event = 1:nr_events,
    
    event_pattern(event_order(event), event:end) = 1;
    
end

events = event_pattern == 1;
events_nonclin = events;
noevents = event_pattern == 0;
noevents_nonclin = noevents;
log_likelihood = zeros(nr_phase, nr_pat);
for phase = 1:nr_phase,
    
    logPrMuEvent=log(prior_thisMean_event(events_nonclin(:, phase)'));
    logPrCovEvent=log(prior_thisCov_event(events_nonclin(:, phase)'));
    
    logPrMuNoEvent=log(prior_thisMean_noevent(noevents_nonclin(:, phase)'));
    logPrCovNoEvent=log(prior_thisCov_noevent(noevents_nonclin(:, phase)'));
    
    logPriorMu(phase,:)=[logPrMuEvent logPrMuNoEvent];
    logPriorCov(phase,:)=[logPrCovEvent logPrCovNoEvent];
    
    logPrEvent=log_likelihood_nonclinevents(events_nonclin(:, phase), :, 2);
    
    logPrNoEvent=log_likelihood_nonclinevents(noevents_nonclin(:, phase), :, 1);
    
    log_likelihood(phase, :)=sum(logPrEvent,1)+sum(logPrNoEvent,1);
    
end



if flag_sum == 2,
    
    log_likelihood_tot = sum(log(sum(exp(log_likelihood), 1))) + sum(logPriorMu(:))+sum(logPriorCov(:));
    
end

if nargout == 2,
    
    class_pat = zeros(nr_pat, 1);
    for pat = 1:nr_pat,
        
        I_pat = find(log_likelihood(:, pat) == max(log_likelihood(:, pat)));
        class_pat(pat) = I_pat(1);
        
    end
    
end

