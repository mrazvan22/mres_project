%% Set up simulated data
clear all
close all

nr_events = 3;
idx_nonclinevent = [1 2 3];
idx_clinevent = [];
event_order(:, 1) = [1 2 3];
event_order(:, 2) = [2 1 3];
fraction_order_event = [0.5 0.5];
event_pattern = zeros(nr_events, nr_events, 2);
for mixt = 1:2,
    
    for event = 1:nr_events,

        event_pattern(event_order(:, mixt) == event, event:end, mixt) = 1;

    end
    
end
nr_controls = 100;
nr_patients = 200;
position_patient = ceil(nr_events*rand(nr_patients, 1));
event_pattern_patient = zeros(nr_events, nr_patients);
for pat = 1:nr_patients,
    
    if pat < fraction_order_event(1)*nr_patients,

        event_pattern_patient(:, pat) = ...
            event_pattern(:, position_patient(pat), 1);
        
    else
        
         event_pattern_patient(:, pat) = ...
            event_pattern(:, position_patient(pat), 2);
        
    end      
    
end

gmix_controls = gmm(1, 1, 'spherical');
gmix_controls.centres = 0.95;
gmix_controls.covars = 0.0001;
gmix_patients = gmm(1, 1, 'spherical');
gmix_patients.centres = 0.8;
gmix_patients.covars = 0.0005;

data_controls = zeros(nr_controls, nr_events);
for event = 1:length(idx_nonclinevent),
    
    data_controls(:, idx_nonclinevent(event)) = gmmsamp(gmix_controls, nr_controls);
    
end

data_patients = zeros(nr_patients, nr_events);
for pat = 1:nr_patients,
    
    for event = 1:length(idx_nonclinevent),
        
        if event_pattern_patient(idx_nonclinevent(event), pat) == 1,
            
            data_patients(pat, idx_nonclinevent(event)) = gmmsamp(gmix_patients, 1);
            
        elseif event_pattern_patient(idx_nonclinevent(event), pat) == 0,
            
            data_patients(pat, idx_nonclinevent(event)) = gmmsamp(gmix_controls, 1);
            
        end
        
    end
    for event = 1:length(idx_clinevent),
            
        data_patients(pat, idx_clinevent(event)) = ...
            event_pattern_patient(idx_clinevent(event), pat);
        
    end
    
end
data_tot = [data_controls; data_patients];
data_range = max(data_tot) - min(data_tot);

%% Compute data Likelihood given event/no event
version_likelihood = 6;
[likelihood, gmix_struct] = ...
    EBDPComputeLikelihood(data_patients, data_controls, version_likelihood);

%% Performing MCMC
parm_mcmc.nr_gradient_ascent = 5;
parm_mcmc.nr_it_gradient_ascent = 1e3;
parm_mcmc.nr_it_burnin = 1e4;
parm_mcmc.nr_it_mcmc = 1e4;
parm_mcmc.interval_display = 1e2;
parm_mcmc.nr_it_check_p_false = 1e2;
parm_mcmc.data_range = data_range;
parm_mcmc.range_p_false = [1e-4 0.2
    1e-4 0.2];
parm_mcmc.std_p_false = [1e-4 1e-4]';

[parm_struct] = EBDPMCMCTest(likelihood_events, parm_mcmc, idx_clinevent);

subplot(131), plot(parm_struct.log_likelihood_mcmc)
subplot(132), imagesc(parm_struct.event_order_mcmc)
subplot(133), plot(parm_struct.p_false_mcmc')

