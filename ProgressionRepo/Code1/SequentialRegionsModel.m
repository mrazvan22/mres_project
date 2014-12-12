
%% Load in some appropriate data.
load fAD_34_LH.mat
%atrophy = data_lh_danny';
load fAD_34_RH.mat
%atrophy = data_rh_danny';
atrophy = [data_lh_danny; data_rh_danny]';

% Load from complete patient data set
load data_jacob_patient.mat
% Generate regional data excluding the hippocampus
atrophy = data_jacc_patient([(end-68):(end-35) (end-33):end], :)';
% Generate regional data including the hippocampus
atrophy = data_jacc_patient([(end-69):end], :)';

% Load from complete control data set
load data_jacob_control.mat
% Generate regional data excluding the hippocampus
atrophy = data_jacc_control([(end-68):(end-35) (end-33):end], :)';
% Generate regional data including the hippocampus
atrophy = data_jacc_control([(end-69):end], :)';

load data_jacob_patient_checked.mat
pjac = data_jacc_patient([(end-69):end], :);
load data_jacob_control_checked.mat
cjac = data_jacc_control([(end-69):end], :);

%load PatientJacManHippo.mat
load PatientJacManHippoChecked.mat
% Left hemisphere only
atrophy = pjac(1:35,:)';
% Right hemisphere only
atrophy = pjac(36:70,:)';
% Both hemispheres
atrophy = pjac';

%load ControlJacManHippo.mat
load ControlJacManHippoChecked.mat
% Left hemisphere only
atrophy = cjac(1:35,:)';
% Right hemisphere only
atrophy = cjac(36:70,:)';
% Both hemispheres
atrophy = cjac';



%% Compare atrophies in two groups.
figure; 
mat = mean(cjac,2);
sat = std(cjac');
errorbar(mat, sat)
hold on;
for i=1:70
scatter(ones(29,1)*i, cjac(i,:), 'b')
end
for i=1:70
scatter(ones(32,1)*i, pjac(i,:), 'r')
end
patmeanat = mean(pjac,2);
patjdiff = patmeanat - mat;


%% Histogram of jacobians in each group
hist(pjac(:), 100, 'FaceColor', 'g')
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w')
hold on;
hist(cjac(:), 100)


%% Create the precedence probability matrix.
[numMeas numEvents] = size(atrophy);
%p = PrecedenceMatrix1(atrophy);
%p = PrecedenceMatrix1r(atrophy);
p = PrecedenceMatrix2(cjac', pjac');



%% Initialize the model

% Construct starting path.
[junk samplePath] = sort(rand(1,numEvents));

% Evaluate likelihood
pathlik = 0;
for i=1:(numEvents-1)
    for j=(i+1):numEvents
        pathlik = pathlik + log(p(samplePath(i), samplePath(j)));
    end
end



%% Run the MCMC
curPath = samplePath;
curlik = pathlik;

burnin = 5000;
interval = 50;
numSamples = 2000;
iterations = interval*numSamples+burnin;
samples = zeros(numSamples, numEvents);
liks = zeros(numSamples, 1);

% Run MCMC
sno = 1;
for it=1:iterations
    
    % Perturb current path
    ind1 = fix(rand(1,1)*numEvents)+1;
    ind2 = ind1;
    while(ind2==ind1)
        ind2 = fix(rand(1,1)*numEvents)+1;
    end
    
    % Switch entries
    newPath = curPath;
    newPath(ind1) = curPath(ind2);
    newPath(ind2) = curPath(ind1);
    
    % Evaluate likelihood
    newlik = 0;
    for i=1:(numEvents-1)
        for j=(i+1):numEvents
            newlik = newlik + log(p(newPath(i), newPath(j)));
        end
    end

    if(newlik>curlik)
        acc = 1;
    else
        a = exp(newlik - curlik);
        acc = rand<a;
    end
    
    if(acc)
        curlik = newlik;
        curPath = newPath;
    end

    if(it>burnin && mod(it-burnin, interval)==0)
        samples(sno,:) = curPath;
        liks(sno) = curlik;
        sno = sno+1;
        display(sprintf('It: %i. Lik: %f', it, curlik));
    end
end


   

%% Test likelihoods of optimal paths for one particular simulation
% Usually ignore this part...

% % optPath(1,:) = [4 6 9 2 8 7 10 1 3 5];
% % optPath(2,:) = [4 9 6 2 8 7 10 1 3 5];
% % optPath(3,:) = [4 6 9 2 8 10 7 1 3 5];
% % optPath(4,:) = [4 9 6 2 8 10 7 1 3 5];
% % optPath(5,:) = [4 6 9 2 8 7 10 3 1 5];
% % optPath(6,:) = [4 9 6 2 8 7 10 3 1 5];
% % optPath(7,:) = [4 6 9 2 8 10 7 3 1 5];
% % optPath(8,:) = [4 9 6 2 8 10 7 3 1 5];
% optPath(1,:) = [1 2 3 4 5];
% optPath(2,:) = [1 4 3 2 5];
% optlik = zeros(1,length(optPath(:,1)));
% for n=1:length(optPath(:,1))
%     for i=1:(numEvents-1)
%         for j=(i+1):numEvents
%             optlik(n) = optlik(n) + log(p(optPath(n,i), optPath(n,j)));
%         end
%     end
% end



%% Look at the likelihood distribution
figure; hist(liks, 100);
%xlim([max(optliks), max(optliks)-5]);
%xlim([max(liks)-100, max(liks)]);


%% Extract just the unique samples in ascending order of likelihood.
[uliks ulikinds n] = unique(liks);
samples(ulikinds,:)


%% Compute the average position of each event
% Next version is better.
% numSamples = length(samples);
% sumord=zeros(1,numEvents);
% for i=1:numSamples
%     [a b] = sort(samples(i,:));
%     sumord = sumord+b;
% end
% avPos = sumord/numSamples
% 


%% Compute the average position of each event weighted by path probabilities
sumord=zeros(1,numEvents);
sumprobs = 0;
% To avoid zero probabilities from numerical truncation, often need to use
% higher likelihoods but such that the relative weightings are the same.
adjuliks = uliks*(-400)/min(uliks);
for i=1:length(ulikinds)
    [a b] = sort(samples(ulikinds(i),:));
    prob = exp(adjuliks(i));
    sumord = sumord+b*prob;
    sumprobs = sumprobs + prob;
end
avPos = sumord/sumprobs


%% Save the results
% Unchecked segmentations and registrations
%save PatientsBothHemiSeqRegionsOutput.mat samples burnin interval numSamples avPos p atrophy
%save PatientsLeftHemiSeqRegionsOutput.mat samples burnin interval numSamples avPos p atrophy
%save PatientsRightHemiSeqRegionsOutput.mat samples burnin interval numSamples avPos p atrophy
%save ControlsBothHemiSeqRegionsOutput.mat samples burnin interval numSamples avPos p atrophy
%save ControlsLeftHemiSeqRegionsOutput.mat samples burnin interval numSamples avPos p atrophy
%save ControlsRightHemiSeqRegionsOutput.mat samples burnin interval numSamples avPos p atrophy

%save PatientsBothHemiSeqRegionsP2_OrigOutput.mat samples burnin interval numSamples avPos p PatAtrophy ConAtrophy
%save PatientsBothHemiSeqRegionsP2_EqualOutput.mat samples burnin interval numSamples avPos p PatAtrophy ConAtrophy patjdiff

% Checked segmentations and registrations
%save PatientsCheckedBothHemiSeqRegionsP1_Output.mat samples burnin interval numSamples avPos p atrophy
%save PatientsCheckedLeftHemiSeqRegionsP1_Output.mat samples burnin interval numSamples avPos p atrophy
%save PatientsCheckedRightHemiSeqRegionsP1_Output.mat samples burnin interval numSamples avPos p atrophy
save ControlsCheckedBothHemiSeqRegionsP1_Output.mat samples burnin interval numSamples avPos p atrophy
%save ControlsCheckedLeftHemiSeqRegionsP1_Output.mat samples burnin interval numSamples avPos p atrophy
%save ControlsCheckedRightHemiSeqRegionsP1_Output.mat samples burnin interval numSamples avPos p atrophy

%save PatientsCheckedBothHemiSeqRegionsP2_OrigOutput.mat samples burnin interval numSamples avPos p PatAtrophy ConAtrophy
%save PatientsCheckedBothHemiSeqRegionsP2_EqualOutput.mat samples burnin interval numSamples avPos p PatAtrophy ConAtrophy patjdiff


%% Plot the average positions with labels

figure; 

% Horizontal positions of each marker are just their average position in
% the order.
% Vertical positions are displaced to avoid overlap of markers and
% emphasize the clustering.  If a marker does not overlap with its
% predecessor, place it at vertical position of zero.  Otherwise move it up
% far enough not to overlap.  For sequences of overlapping ones, we move
% them alternately up and down.
[ord inds] = sort(avPos);
vertCo = zeros(size(avPos));
% For one hemisphere
%xRad = 0.5;
%yRad = 0.04;
% For both hemispheres
xRad = 1.2;
yRad = 0.3;
for i=2:length(avPos)
    dist = avPos(inds(i))-avPos(inds(i-1));
    if(dist<2*xRad)
        %disp = 2*yRad*sqrt(1-(dist/(4*xRad^2)));
        if(vertCo(inds(i-1))>=0)
            vertCo(inds(i)) = -(vertCo(inds(i-1)) + 2*yRad);
        else
            vertCo(inds(i)) = -vertCo(inds(i-1));
        end
    end
end

for i=1:35
    labels{i} = sprintf('%02iL', i);
    labels{i+35} = sprintf('%02iR', i);
end
% One hemisphere
%scatter(avPos, vertCo, 'SizeData', 500, 'LineWidth', 3);
% Both hemispheres
%scatter(avPos, vertCo, 'SizeData', 1000, 'LineWidth', 3);
% This is for two-sided probability and distinguishes expanders from
% shrinkers.
%Left
%scatter(avPos, vertCo, 500*ones(length(avPos),1), patjdiff(1:35), 'SizeData', 500, 'LineWidth', 3);
%Right
%scatter(avPos, vertCo, 500*ones(length(avPos),1), patjdiff(36:70), 'SizeData', 500, 'LineWidth', 3);
%Both
scatter(avPos, vertCo, 1000*ones(length(avPos),1), patjdiff, 'SizeData', 1000, 'LineWidth', 3);
colormap winter;
hold on;
for i=1:length(avPos);
    % This one is for both hemispheres (68 or 70 regions)
    text(avPos(i)-0.7, vertCo(i), labels{i});
    % This one works for the 34 or 35 regions (single hemisphere)
    %text(avPos(i)-0.2, vertCo(i), labels{i}(1:2));
end
% One hemi
%ylim([min(vertCo),-min(vertCo)]*5/2);
% Both hemis
ylim([min(vertCo),-min(vertCo)+0.2]);

axis off;

%% Write out the figure
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 14 6]);
%print -dpng 'PatientsBothHemiSeqRegionsOrdering.png'
%print -dpng 'PatientsLeftHemiSeqRegionsOrdering.png'
%print -dpng 'PatientsRightHemiSeqRegionsOrdering.png'
%print -dpng 'ControlsBothHemiSeqRegionsOrdering.png'
%print -dpng 'ControlsLeftHemiSeqRegionsOrdering.png'
%print -dpng 'ControlsRightHemiSeqRegionsOrdering.png'
%print -dpng 'PatientsBothHemiSeqRegionsP2_OrigOrdering.png'
%print -dpng 'PatientsBothHemiSeqRegionsP2_EqualOrdering.png'

%print -dpng 'PatientsCheckedBothHemiSeqRegionsP1_Ordering.png'
%print -dpng 'PatientsCheckedLeftHemiSeqRegionsP1_Ordering.png'
%print -dpng 'PatientsCheckedRightHemiSeqRegionsP1_Ordering.png'
print -dpng 'ControlsCheckedBothHemiSeqRegionsP1_Ordering.png'
%print -dpng 'ControlsCheckedLeftHemiSeqRegionsP1_Ordering.png'
%print -dpng 'ControlsCheckedRightHemiSeqRegionsP1_Ordering.png'
%print -dpng 'PatientsCheckedBothHemiSeqRegionsP2_OrigOrdering.png'
%print -dpng 'PatientsCheckedBothHemiSeqRegionsP2_EqualOrdering.png'


%% Reorder the list of region labels
load roi_name_select.mat
RegionLabels{1} = 'ctx-lh-hippocampus';
RegionLabels{36} = 'ctx-rh-hippocampus';
for i=2:35
    RegionLabels{i} = label_select{i+1};
    RegionLabels{i+35} = label_select{i+35};
end
save RegionLabels.mat RegionLabels;


%% Make a figure of the key
figure;
load RegionLabels.mat
rows = 6;
for i=1:(length(label_select)/2)
    text(fix((i-1)/rows)*200, rows-1-mod(i-1,rows), sprintf('%02i\t%s', i, RegionLabels{i}(8:end)));
end
xlim([0 200*(rows+1)]);
ylim([0, rows+1]);
axis off;
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 24 3]);
print -dpng 'LabelKey.png'



%% Alternative very simple ordering from average atrophy in each region
avPos = mean(atrophy);
avPos = (avPos-min(avPos))/(max(avPos)-min(avPos))*length(avPos);
% This produces a similar ordering to the above.

%% Write out the figure
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 14 6]);
print -dpng 'PatientsCheckedBothHemiMeanAtrophyOrdering.png'



%% Compare the orderings in patients and controls
load ControlsBothHemiSeqRegionsOutput.mat
AvPosControl = avPos;
load PatientsBothHemiSeqRegionsOutput.mat
AvPosPatient = avPos;

figure; 
scatter(AvPosPatient, AvPosControl, 'SizeData', 500, 'LineWidth', 3);
hold on;
for i=1:length(AvPosControl)
    text(AvPosPatient(i)-2, AvPosControl(i), labels{i});
end
xlabel('Patient av. pos.')
ylabel('Control av. pos.')

%% Write out the figure
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6 6]);
print -dpng 'PatientVControlsBothHemiSeqRegionsOrdering.png'


%% Compare orderings in patient left and right hemispheres.
figure; 
scatter(AvPosPatient(1:35), AvPosPatient(36:70), 'SizeData', 500, 'LineWidth', 3);
hold on;
for i=1:(length(AvPosPatient)/2)
    text(AvPosPatient(i)-2, AvPosPatient(i+35), labels{i});
end
xlabel('Patient Left av. pos.')
ylabel('Patient Right av. pos.')

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6 6]);
print -dpng 'PatientLeftVRightSeqRegionsOrdering.png'

%% Now for controls
figure; 
scatter(AvPosControl(1:35), AvPosControl(36:70), 'SizeData', 500, 'LineWidth', 3);
hold on;
for i=1:(length(AvPosControl)/2)
    text(AvPosControl(i)-2, AvPosControl(i+35), labels{i});
end
xlabel('Control Left av. pos.')
ylabel('Control Right av. pos.')

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6 6]);
print -dpng 'ControlLeftVRightSeqRegionsOrdering.png'


%% Compute precedence matrix based on comparison of patients with control distribution.
%load data_jacob_control.mat
%ConAtrophy = data_jacc_control([(end-69):end], :)';
load ControlJacManHippo.mat
ConAtrophy = cjac';
mat = mean(ConAtrophy);
sat = std(ConAtrophy);
mmat = mean(ConAtrophy(:));
msat = std(ConAtrophy(:));

%load data_jacob_patient.mat
%PatAtrophy = data_jacc_patient([(end-69):end], :)';
load PatientJacManHippo.mat
PatAtrophy = pjac';
[numMeas numEvents] = size(PatAtrophy);

% Compare atrophies in two groups.
figure; 
errorbar(mat, sat)
hold on;
for i=1:70
scatter(ones(29,1)*i, ConAtrophy(:,i), 'b')
end
for i=1:70
scatter(ones(32,1)*i, PatAtrophy(:,i), 'r')
end
patmeanat = mean(PatAtrophy);
patjdiff = patmeanat - mat;

% Probability of atrophy for each region in each patient data set assuming
% uniform distribution for jacobian given atrophy and using the normal
% model of unatrophied jacobians from controls.
pat = zeros(size(PatAtrophy));
for i=1:numMeas
    for j=1:numEvents
        % Regional means and standard deviations
%        tp = exp(-(PatAtrophy(i,j) - mat(j))^2/(2*sat(j)^2))/(sqrt(2*pi*sat(j)^2));
        % Pooled mean and standard deviation
        tp = exp(-(PatAtrophy(i,j) - mmat)^2/(2*msat^2))/(sqrt(2*pi*msat^2));
        pat(i,j) = 1/(1+tp);
    end
end

% Probability of atrophy in the first of a pair, but not the second in each
% patient data set.
par1nar2 = zeros(numMeas, numEvents, numEvents);
for i=1:numMeas
    par1nar2(i,:,:) = pat(i,:)'*(1-pat(i,:));
end


% Now produce the p(R1<R2) matrix.
% May need to be careful about entries for which probability is close to
% zero everywhere, but this does not seem to be a problem at the moment.
% If so, add a regularizing factor to the denominator.
p = squeeze(sum(par1nar2.^2)./sum(par1nar2));

% Probabilities of simultaneous atrophy.  Here we use a simple but naive
% definition, but it should really account better for the distribution of
% samples along the L-curve.
%peq = (1 - p).*(1 - p');

