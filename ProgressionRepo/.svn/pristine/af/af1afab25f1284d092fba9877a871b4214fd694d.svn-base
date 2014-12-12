%% Simulation of AD atrophy in different regions

%% Quadratic regional atrophy starting at different times in different regions.
quad_hippovol = @(t)(1-t.^2);
quad_precuvol = @(t)(1-(max([0.2*ones(length(t),1),t']')-0.2).^2);
quad_tlcvol = @(t)(1-(max([0.4*ones(length(t),1),t']')-0.4).^2);
quad_pcgvol = @(t)(1-(max([0.4*ones(length(t),1),t']')-0.4).^2);
quad_accvol = @(t)(1-(max([0.7*ones(length(t),1),t']')-0.7).^2);

% Plot them for comparison.
figure;
t=0:0.01:1;
plot(t, quad_hippovol(t));
hold on;
plot(t, quad_precuvol(t),'g');
plot(t, quad_tlcvol(t),'r');
plot(t, quad_pcgvol(t),'m');
plot(t, quad_accvol(t),'c');


%% Linear regional atrophy starting at different times in different regions.
lin_hippovol = @(t)(1-t);
lin_precuvol = @(t)(1-(max([0.2*ones(length(t),1),t']')-0.2));
lin_tlcvol = @(t)(1-(max([0.4*ones(length(t),1),t']')-0.4));
lin_pcgvol = @(t)(1-(max([0.4*ones(length(t),1),t']')-0.4));
lin_accvol = @(t)(1-(max([0.7*ones(length(t),1),t']')-0.7));

% Plot them for comparison.
figure;
t=0:0.01:1;
plot(t, lin_hippovol(t));
hold on;
plot(t, lin_precuvol(t),'g');
plot(t, lin_tlcvol(t),'r');
plot(t, lin_pcgvol(t),'m');
plot(t, lin_accvol(t),'c');


%% Get measurements of relative volume at lots of random time points.
NumSubjects = 100;

stimes = rand(1,NumSubjects);

%rvols = [lin_hippovol(stimes); lin_precuvol(stimes); lin_tlcvol(stimes); lin_pcgvol(stimes); lin_accvol(stimes)]';
rvols = [quad_hippovol(stimes); quad_precuvol(stimes); quad_tlcvol(stimes); quad_pcgvol(stimes); quad_accvol(stimes)]';
[numMeas numEvents] = size(rvols);

% Add a little bit of noise
rvols = rvols+randn(numMeas, numEvents)*0.01;


%% Write out a file in the format required by biolearn.
%fid = fopen('AD_Sim.txt', 'w');
%fprintf(fid, 'HIPP\tPREC\tTLC\tPCG\tACC\n');
%for i=1:numMeas
%    fprintf(fid, '%d\t%d\t%d\t%d\t%d\n', rvols(i,1), rvols(i,2), rvols(i,3), rvols(i,4), rvols(i,5));
%end
%fclose(fid);


%% Simulate jacobians from fAD data by drawing both baseline and follow up scan times randomly along the progression axis.

% Assume baseline images are acquired in presymtomatic phase, which we'll
% say is between 0 and 0.4
baselines = rand(1,NumSubjects)*0.4;

% Follow ups occur at random times between baseline and the end of the
% disease.
followups = rand(1,NumSubjects).*(1-baselines)+baselines;

% Compute volumes for each region at each timepoint.
blvols = [hippovol(baselines); precuvol(baselines); tlcvol(baselines); pcgvol(baselines); accvol(baselines)];
fuvols = [hippovol(followups); precuvol(followups); tlcvol(followups); pcgvol(followups); accvol(followups)];
NumRegions = 5;

% Compute regional jacobians
pjac = fuvols./blvols;

% Create artificial control jacobians
cjac = randn(size(pjac))*0.05 + 1;


%% Compare atrophies in two groups.
figure; 
mat = mean(cjac,2);
sat = std(cjac');
errorbar(mat, sat)
hold on;
for i=1:NumRegions
scatter(ones(NumSubjects,1)*i, cjac(i,:), 'b')
end
for i=1:NumRegions
scatter(ones(NumSubjects,1)*i, pjac(i,:), 'r')
end
patmeanat = mean(pjac,2);
patjdiff = patmeanat - mat;


%% Estimate the precedence matrix.

p = PrecedenceMatrixZ2(cjac', pjac');
figure; imagesc(p);


