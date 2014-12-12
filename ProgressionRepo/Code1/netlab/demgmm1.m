%DEMGMM1 Demonstrate EM for Gaussian mixtures.
%
%	Description
%	This script demonstrates the use of the EM algorithm to fit a mixture
%	of Gaussians to a set of data using maximum likelihood. A colour
%	coding scheme is used to illustrate the evaluation of the posterior
%	probabilities in the E-step of the EM algorithm.
%
%	See also
%	DEMGMM2, DEMGMM3, DEMGMM4, GMM, GMMEM, GMMPOST
%

%	Copyright (c) Ian T Nabney (1996-2001)

clc;
disp('This demonstration illustrates the use of the EM (expectation-')
disp('maximization) algorithm for fitting of a mixture of Gaussians to a')
disp('data set by maximum likelihood.')
disp(' ')
disp('The data set consists of 40 data points in a 2-dimensional')
disp('space, generated by sampling from a mixture of 2 Gaussian')
disp('distributions.')
disp(' ')
disp('Press any key to see a plot of the data.')
pause;

% Generate the data
randn('state', 0); rand('state', 0);
gmix = gmm(2, 2, 'spherical');
ndat1 = 20; ndat2 = 20; ndata = ndat1+ndat2;
gmix.centres =  [0.75 0.75; 0.25 0.25]; 
gmix.covars = [0.01 0.01];
x = gmmsamp(gmix, ndata);

h = figure;
hd = plot(x(:, 1), x(:, 2), '.g', 'markersize', 30);
hold on; axis([0 1 0 1]); axis square; set(gca, 'box', 'on');
ht = text(0.5, 1.05, 'Data', 'horizontalalignment', 'center');
disp(' ');
disp('Press any key to continue.')
pause; clc;

disp('We next create and initialize a mixture model consisting of a mixture')
disp('of 2 Gaussians having ''spherical'' covariance matrices, using the')
disp('function GMM. The Gaussian components can be displayed on the same')
disp('plot as the data by drawing a contour of constant probability density')
disp('for each component having radius equal to the corresponding standard')
disp('deviation. Component 1 is coloured red and component 2 is coloured')
disp('blue.')
disp(' ')
disp('Note that a particulary poor choice of initial parameters has been')
disp('made in order to illustrate more effectively the operation of the')
disp('EM algorithm.')
disp(' ')
disp('Press any key to see the initial configuration of the mixture model.')
pause;

% Set up mixture model
ncentres = 2; input_dim = 2;
mix = gmm(input_dim, ncentres, 'spherical');

% Initialise the mixture model
mix.centres = [0.2 0.8; 0.8, 0.2];
mix.covars = [0.01 0.01];

% Plot the initial model
ncirc = 30; theta = linspace(0, 2*pi, ncirc);
xs = cos(theta); ys = sin(theta);
xvals = mix.centres(:, 1)*ones(1,ncirc) + sqrt(mix.covars')*xs;
yvals = mix.centres(:, 2)*ones(1,ncirc) + sqrt(mix.covars')*ys;
hc(1)=line(xvals(1,:), yvals(1,:), 'color', 'r');
hc(2)=line(xvals(2,:), yvals(2,:), 'color', 'b');
set(ht, 'string', 'Initial Configuration');
figure(h);
disp(' ')
disp('Press any key to continue'); 
pause; clc;

disp('Now we adapt the parameters of the mixture model iteratively using the')
disp('EM algorithm. Each cycle of the EM algorithm consists of an E-step')
disp('followed by an M-step.  We start with the E-step, which involves the')
disp('evaluation of the posterior probabilities (responsibilities) which the')
disp('two components have for each of the data points.')
disp(' ')
disp('Since we have labelled the two components using the colours red and')
disp('blue, a convenient way to indicate the value of a posterior')
disp('probability for a given data point is to colour the point using a')
disp('scale ranging from pure red (corresponding to a posterior probability')
disp('of 1.0 for the red component and 0.0 for the blue component) through')
disp('to pure blue.')
disp(' ')
disp('Press any key to see the result of applying the first E-step.')
pause;

% Initial E-step.
set(ht, 'string', 'E-step');
post = gmmpost(mix, x);
dcols = [post(:,1), zeros(ndata, 1), post(:,2)];
delete(hd); 
for i = 1 : ndata
  hd(i) = plot(x(i, 1), x(i, 2), 'color', dcols(i,:), ...
          'marker', '.', 'markersize', 30);
end
figure(h);

disp(' ');
disp('Press any key to continue')
pause; clc;

disp('Next we perform the corresponding M-step. This involves replacing the')
disp('centres of the component Gaussians by the corresponding weighted means')
disp('of the data. Thus the centre of the red component is replaced by the')
disp('mean of the data set, in which each data point is weighted according to')
disp('the amount of red ink (corresponding to the responsibility of')
disp('component 1 for explaining that data point). The variances and mixing')
disp('proportions of the two components are similarly re-estimated.')
disp(' ')
disp('Press any key to see the result of applying the first M-step.')
pause;

% M-step.
set(ht, 'string', 'M-step');
options = foptions; 
options(14) = 1; % A single iteration
options(1) = -1; % Switch off all messages, including warning
mix = gmmem(mix, x, options);
delete(hc);
xvals = mix.centres(:, 1)*ones(1,ncirc) + sqrt(mix.covars')*xs;
yvals = mix.centres(:, 2)*ones(1,ncirc) + sqrt(mix.covars')*ys;
hc(1)=line(xvals(1,:), yvals(1,:), 'color', 'r');
hc(2)=line(xvals(2,:), yvals(2,:), 'color', 'b');
figure(h);
disp(' ')
disp('Press any key to continue')
pause; clc;

disp('We can continue making alternate E and M steps until the changes in')
disp('the log likelihood at each cycle become sufficiently small.')
disp(' ')
disp('Press any key to see an animation of a further 9 EM cycles.')
pause;
figure(h);

% Loop over EM iterations.
numiters = 9;
for n = 1 : numiters

  set(ht, 'string', 'E-step');
  post = gmmpost(mix, x);
  disp(post)
  dcols = [post(:,1), zeros(ndata, 1), post(:,2)];
  delete(hd); 
  for i = 1 : ndata
    hd(i) = plot(x(i, 1), x(i, 2), 'color', dcols(i,:), ...
                 'marker', '.', 'markersize', 30);
  end
  pause(1)

  set(ht, 'string', 'M-step');
  [mix, options] = gmmem(mix, x, options);
  fprintf(1, 'Cycle %4d  Error %11.6f\n', n, options(8));
  delete(hc);
  xvals = mix.centres(:, 1)*ones(1,ncirc) + sqrt(mix.covars')*xs;
  yvals = mix.centres(:, 2)*ones(1,ncirc) + sqrt(mix.covars')*ys;
  hc(1)=line(xvals(1,:), yvals(1,:), 'color', 'r');
  hc(2)=line(xvals(2,:), yvals(2,:), 'color', 'b');
  pause(1)

end

disp(' ');
disp('Press any key to end.')
pause; clc; close(h); clear all

