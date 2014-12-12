%% Making mixture illustration

gmix_controls = gmm(1, 1, 'spherical');
gmix_controls.centres = 0.95;
gmix_controls.covars = 0.002;
gmix_patients = gmm(1, 2, 'spherical');
gmix_patients.centres = [0.95 0.75]';
gmix_patients.covars = [0.002 0.002];
gmix_patients.priors = [0.2 0.8];
gmix_outliers = gmm(1, 1, 'spherical');
gmix_outliers.centres = 0.75;
gmix_outliers.covars = [0.002];

data_controls = gmmsamp(gmix_controls, 200);
data_patients = gmmsamp(gmix_patients, 200);
[n_c, x_c] = hist(data_controls, 10);
[n_p, x_p] = hist(data_patients, 20);
n_c_max = max(n_c);
n_p_max = max(n_p);

x_min = min(min(data_controls), min(data_patients));
x_max = max(max(data_controls), max(data_patients));
x_prob = linspace(x_min, x_max, 1e3)';
Lik_controls = gmmprob(gmix_controls, x_prob);
Lik_outliers = gmmprob(gmix_outliers, x_prob);
Lik_controls = Lik_controls*((1.2*n_c_max)/max(Lik_controls));
Lik_outliers = Lik_outliers*((1.2*n_p_max)/max(Lik_outliers));

figure(1), clf, hold on
bar(x_c, n_c, 'g');
bar(x_p, n_p, 'r');
ax1 = gca;
hxaxis1 = xlabel('event measure (a.u.)');
hyaxis1 = ylabel('number of subjects');
ax2 = axes('Position', get(ax1,'Position'),...
           'XAxisLocation','bottom',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
hl21 = line(x_prob, Lik_controls, 'Color', 'g', 'LineWidth', 4, 'Parent', ax2);
hl22 = line(x_prob, Lik_outliers, 'Color', 'r', 'LineWidth', 4, 'Parent', ax2);
hyaxis2 = ylabel('Likelihood');
set([hxaxis1, hyaxis1, hyaxis2], ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

%% Likelihood figures

clear all
close all

event = rand(4, 8);
event = event > 0.5;



