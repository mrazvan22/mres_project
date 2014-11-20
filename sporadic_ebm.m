function [] = sporadic_ebm()

% ADNIdataBL raw data as received from the ADNI website
%EBMdataBL - pre-processed data for use in the EBM model
% EMBLdxBL - the level of impairment for each patient (1 - Cognitively normal, 2 - MCI 3 - AD)
% EBMevents - labels of the EBM events
load('alex_data/ADNIdata_Baseline.mat')

%[mu_mix, sigma_mix, pi_mix] = calc_gaussian_parameters(EBMdataBL, EBMdxBL)
load('sporadic_params.mat')

%grad_asc_S = char_seq_grad_asc(EBMdataBL, mu_mix, sigma_mix)
%save('sporadic_params.mat', 'grad_asc_S','-append')
MCMC_samples = MCMC_sampling(EBMdataBL, mu_mix, sigma_mix, grad_asc_S);

calc_variance_diagrams(samples);

%calc_likelihood(EBMdataBL, 1:14, mu_mix, sigma_mix)

end

function [mu_mix, sigma_mix, pi_mix] = calc_gaussian_parameters(EBMdataBL, EBMdxBL)


[nr_patients, nr_biomarkers] = size(EBMdataBL);

control_indices = find(EBMdxBL == 1);
patient_indices = find(EBMdxBL > 1);

% pXgEnE(i,j,1) = p(x_ij | Ei) (patients)
% pXgEnE(i,j,2) = p(x_ij | not Ei) (controls)
pXgEnE = zeros(nr_patients, nr_biomarkers,2);

mu_mix = zeros(nr_biomarkers, 2);
sigma_mix = zeros(nr_biomarkers, 2);
pi_mix = zeros(nr_biomarkers, 2);

for biomk=1:nr_biomarkers
   control_levels = EBMdataBL(control_indices, biomk);
   patient_levels = EBMdataBL(patient_indices, biomk);
    
   mu_control = mean(control_levels);
   sigma_control = std(control_levels);
   
   mu_patient = mean(patient_levels);
   sigma_patient = std(patient_levels);
   
   %pXgEnE(:, biomk, 2) = normpdf(control_levels, mu_control, sigma_control);
   %pXgEnE(:, biomk, 1) = normpdf(patient_levels, mu_patient, sigma_patient);
   
   
   all_levels = EBMdataBL(:, biomk);
   % fit two gaussians on all the data

   [mu_mix(biomk,:), sigma_mix(biomk,:), pi_mix(biomk,:)]  = ...
       em_mix(all_levels, [mu_control, mu_patient], [sigma_control, sigma_patient])
   

   minX = min(all_levels);
   maxX = max(all_levels);
   X = minX:(maxX - minX)/200:maxX;
   Y = eval_mix_gaussians(X, mu_mix(biomk,:), sigma_mix(biomk,:), pi_mix(biomk,:));
   %Y = eval_mix_gaussians(X, [mu_control, mu_patient], [sigma_control, sigma_patient], pi_mix(biomk,:));
   
   
   clf;
   [f1,x1] = hist(patient_levels);
   bar(x1,f1, 'FaceColor',[0.8,0.2,0.2]);
   %h = findobj(gca,'Type','patch');
   %h.FaceColor = [0.8 0.2 0.2];
   %h.EdgeColor = 'b';
   hold on
   [f2,x2] = hist(control_levels);
   bar(x2,f2,'b');
   maxY = get(gca,'ylim');
   Y = Y * (maxY(2) * 0.66/max(Y));
   %hist(control_levels);
   hold on

   plot(X, Y,'y');

   %hist(control_levels,50)
   
   if(biomk == 14)
        display('asfda')
   end
   
end

   save('sporadic_params.mat', 'mu_mix', 'sigma_mix', 'pi_mix'); 
   

end

function [mu, sigma, pi]  = ...
    em_mix(data, mu_init, max_sigma)

% data is a 1d array of biomarker levels which is assumed to be generated
% from a mixture of 2 gaussian distributions

mu = mu_init;
sigma = max_sigma;
min_sigma = max_sigma;

% pi_control is the weight of the control gaussian
% pi(1) - control pi(2) - patient
pi = [0.5 0.5];

N = length(data); 
K = 2;

assert(length(mu) == K && length(sigma) == K && length(pi) == K);

znk = zeros(N,K);
eval_matrix = zeros(N,K);

log_likely = inf;
thresh = 0.0001;

iterations = 300;

Nk = [-1 -1];
%% keep updating until convergence
for iter=1:iterations
    
    % E-step
    for k=1:K
        znk(:,k) = pi(k) * normpdf(data, mu(k), sigma(k));
    end
% 
    % normalise the znk rows
    for n=1:N
        if(sum(znk(n,:)) ~= 0)
            znk(n,:) = znk(n,:) ./ sum(znk(n,:));
        else
            display('Warning - znk is zero')
            znk(n,:) = [0.5 0.5];
        end
    end

    % M-step

    Nk = sum(znk);
    
    if(~all(Nk) || isnan(Nk(1)) || isnan(Nk(2)))
        break
    end
    
    for k=1:K
        
       mu(k) =  sum(znk(:,k) .* data) /Nk(k);
       sigma(k) = sum(znk(:,k) .* (data - mu(k)) .* (data - mu(k)) ) / Nk(k);

       % set the constraint that alex implemented in the paper: i.e. the sigma
       % of each distribution should not be greater than the sigma of the CN
       % and AD groups taken separately
       if (sigma(k) > max_sigma(k))
           display('max sigma constraint')
           sigma(k) = max_sigma(k);
       end
      
       % also make sure the sigma doesn't fall below a certain threshold
       if (sigma(k) < min_sigma(k))
           display('min sigma constraint')
           sigma(k) = min_sigma(k);
       end
       
       pi = Nk / N;
       
       % normalise the pi, although it shouldn't normally need to be
       % normalised
       pi = pi ./ sum(pi);

       %X = min(data):0.005:max(data);
       %hist(data,25)
       %hold on
       %plot(X, eval_mix_gaussians(X, mu, sigma, pi))

    end
    
           
    % Evaluate the log_likelihood
    for j=1:K
       eval_matrix(:,j) = pi(j) * normpdf(data, mu(j), sigma(j));
    end

    sumK = sum(eval_matrix, 2);
    assert(length(sumK) == N);
    log_likely = sum(log(sumK));
    
    if(~all(sigma) || isnan(sigma(1)) || isnan(sigma(2)))
        break
    end
end
    

end

function logL = calc_likelihood(X, S, mu_mix, sigma_mix)

% modelled from equation (1), alex paper page 2567
% J - # patients
% I - # Events or disease stages
[J,I] = size(X);


[I2,~] = size(mu_mix);
[I3,~] = size(sigma_mix);

assert(I == I2 && I == I3);

logL = 1;

%reorder the gaussian parameters according to the S ordering provided
mu_mix = mu_mix(S,:);
sigma_mix = sigma_mix(S,:);

%also reorder the dataset
X = X(:,S);

% x_{ij} biomarker i in subject j
% pXgE(i,j,1) = p (x_{ij} | E_i)      patients
% pXgE(i,j,2) = p (x_{ij} | not E_i)  controls
pXgE = zeros(I, J, 2);

% for each biomarker 
% mu_mix(:,1) - controls
% mu_mix(:,2) - patients
for biomk=1:I
    pXgE(biomk,:,1) = normpdf(X(:,biomk), mu_mix(biomk,2), sigma_mix(biomk,2));
    pXgE(biomk,:,2) = normpdf(X(:,biomk), mu_mix(biomk,1), sigma_mix(biomk,1));
    if(~all(pXgE(biomk,:,1)))
        display('pXgE is zero')
    end
end

pK = 1/J; % uniform prior that patient i is at stage k.

% for each patient
for patient=1:J
   sum = 0;
   % for each disease stage
   for stage=0:I
        
        sum = sum + prod(pXgE(1:stage,patient,1)) * prod(pXgE(stage+1:I,patient,2));
   end
   if(sum == 0)
      display('Warning log-likelihood is -inf') 
   end
   
   logL = logL  + log(sum * pK);
   
end



end

function curr_seq = char_seq_grad_asc(X, mu_mix, sigma_mix)

[nr_subjects,nr_biomk] = size(X);


[nr_biomk2,~] = size(mu_mix);
[nr_biomk3,~] = size(sigma_mix);

assert(nr_biomk == nr_biomk2 && nr_biomk == nr_biomk3);

curr_seq = 1:14;
nr_iterations = 2000;
old_likelihood = calc_likelihood(X, curr_seq, mu_mix, sigma_mix);

s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);

for i=1:nr_iterations
    i
    p1 = ceil(rand * 14);
    p2 = ceil(rand * 14);
    while (p1 == p2)
        p2 = ceil(rand * 14);
    end
    new_seq = curr_seq;
    tmp = new_seq(p1);
    new_seq(p1) = new_seq(p2);
    new_seq(p2) = tmp;
    new_likelihood = calc_likelihood(X, new_seq, mu_mix, sigma_mix);
    if(new_likelihood >= old_likelihood)
        curr_seq = new_seq;
        old_likelihood = new_likelihood;
    end
    
    curr_seq
    
end


end

% generic function for
function samples = MCMC_sampling(X, mu_mix, sigma_mix, init_seq)

[nr_subjects,nr_biomk] = size(X);


[nr_biomk2,~] = size(mu_mix);
[nr_biomk3,~] = size(sigma_mix);

assert(nr_biomk == nr_biomk2 && nr_biomk == nr_biomk3);

% initialise curr_seq to initial sequence
curr_seq = init_seq;

%burnout_iterations = 10^5;
burnout_iterations = 100;
actual_iterations = 10^6;

old_likelihood = calc_likelihood(X, curr_seq, mu_mix, sigma_mix);

s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);

samples = zeros(actual_iterations, nr_biomk);

for i=1:(burnout_iterations + actual_iterations)
    i
    p1 = ceil(rand * 14);
    p2 = ceil(rand * 14);
    while (p1 == p2)
        p2 = ceil(rand * 14);
    end
    new_seq = curr_seq;
    tmp = new_seq(p1);
    new_seq(p1) = new_seq(p2);
    new_seq(p2) = tmp;
    new_likelihood = calc_likelihood(X, new_seq, mu_mix, sigma_mix);
    likely_ratio = new_likelihood / old_likelihood;
    
   
    if(rand < likely_ratio)
        % swapt them with probability min(likely_ratio, 1)
        curr_seq = new_seq;
        old_likelihood = new_likelihood;
    end
    
    if (i > burnout_iterations)
        samples(i-burnout_iterations,:) = curr_seq;
    end
    
    if(mod(i,1000) == 0)
       display('1000 iterations') 
    end
    
end


end

% function l = get_likelyhood_ratio(new_likelihood, old_likelihood, mode)
% 
% switch mode
%     
%     case 'MCMC'
%         l = new_likelihood / old_likelihood;
%     case 'greedy'
%         if (new_likelihood > old_likelihood)
%             l = 1;
%         else
%             l = 0;
%         end
%     otherwise
%         error('mode has to be MCMC or greedy')
% end
%     
% end

function Y = eval_mix_gaussians(X, mu, sigma, pi)
% X is a vector of data points
% mu is a vector of means
% sigma is a vector of std deviations

K = length(mu);
Y = zeros(size(X));

for k=1:K
   if (sigma(k) ~= 0) 
       Y = Y + pi(k) * normpdf(X, mu(k), sigma(k)); 
   end
end

end

function [] = calc_variance_diagrams(samples)
    
biomk_avg_ordering = mean(MCMC_samples);
[~, char_sequence] = sort(biomk_avg_ordering);
char_sequence


    
end
