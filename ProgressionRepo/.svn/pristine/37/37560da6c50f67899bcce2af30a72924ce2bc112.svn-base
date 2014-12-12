function VerifyMCMCWithParamsDemo
clear
close all

%% Initialize the Gaussian parameters
rand('seed',14);
randn('seed',14);

gauss(1).mean=0;
gauss(1).cov=1;
gauss(1).prior=0.5;

gauss(2).mean=-1;
gauss(2).cov=1;
gauss(2).prior=0.5;

nData=1000;
prGaussCum = cumsum([gauss(:).prior]);


%% Gerenate data
% for (cData = 1:nData)
%     %choose which Gaussian
%     randVal = rand(1);
%     decision = prGaussCum-randVal;
%     %retain gaussians where cum prob is greater
%     decision(find(decision<0))=1.1;
%     %choose closest Gaussian
%     gaussChosen = find(decision == min(decision));
%     %generate this point and store
%     data(cData) = gauss(gaussChosen).mean+gauss(gaussChosen).cov*randn;
%     membership(cData)=gaussChosen;
% end;
% %% Data Visualization
% data_controls=data(find(membership==1));
% data_patients=data(find(membership==2));

data_controls=gauss(1).mean+gauss(1).cov*randn(1,nData/2);
data_patients=gauss(2).mean+gauss(2).cov*randn(1,nData/2);


figure;
[hist_c, x_c] = ksdensity(data_controls);
plot(x_c, hist_c, 'g');
hold on
[hist_p, x_p] = ksdensity(data_patients);
plot(x_p, hist_p, 'r');
set(gcf,'Color',[1 1 1])

%% Performing MCMC

% Performing mcmc
parm_mcmc.nr_gradient_ascent = 5;
parm_mcmc.nr_it_gradient_ascent = 2e3;
parm_mcmc.nr_it_burnin = 1e4;
parm_mcmc.nr_it_mcmc = 1e4;
parm_mcmc.interval_display = 1e3;
parm_mcmc.flag_sum = 2; % This should always be set to 2.
parm_mcmc.nr_it_check_p_false = 1e2;
parm_mcmc.range_p_false = [0.05 0.1
    0.05 0.1];
parm_mcmc.std_p_false = [1e-4 1e-4]';
parm_mcmc.idx_clinevent = [];


parm_struct = VerifyTestWithParamsFast(data_controls, ... 
    data_patients, parm_mcmc);
 