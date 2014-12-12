
%% Load in some appropriate data.
load PatientJacManHippo.mat
% Left hemisphere only
atrophy = pjac(1:35,:)';
% Right hemisphere only
atrophy = pjac(36:70,:)';
% Both hemispheres
atrophy = pjac';

load ControlJacManHippo.mat
% Left hemisphere only
atrophy = cjac(1:35,:)';
% Right hemisphere only
atrophy = cjac(36:70,:)';
% Both hemispheres
atrophy = cjac';


%% Create the precedence probability matrix.
[numMeas numEvents] = size(atrophy);
p = PrecedenceMatrix1(atrophy);


%% Initialize the model
numGroups = 8;

% Assign the regions randomly to groups
groupLabels = fix(rand(1,numEvents)*numGroups)+1;

% Evaluate likelihood
% Between group likelihood
groupordlik = 0;
logp = log(p);
logmp = log(1-p);
logmp(find(logmp<-15)) = -15;
for i=1:(numEvents-1)
    for j=(i+1):numEvents
        if(groupLabels(i)<groupLabels(j))
            groupordlik = groupordlik + logp(i,j) + logmp(j,i);
        elseif(groupLabels(i)==groupLabels(j))
            groupordlik = groupordlik + logmp(i,j) + logmp(j,i);
        else
            groupordlik = groupordlik + logp(j,i) + logmp(i,j);
        end
    end
end



%% Run the MCMC
curGroups = groupLabels;
curlik = groupordlik;

burnin = 5000;
interval = 50;
numSamples = 5000;
iterations = interval*numSamples+burnin;
samples = zeros(numSamples, numEvents);
liks = zeros(numSamples, 1);

% Run MCMC
sno = 1;
for it=1:iterations
    
    % Perturb current grouping
    ind = fix(rand(1,1)*numEvents)+1;
    grp1 = curGroups(ind);
    grp2 = grp1;
    while(grp2==grp1)
        grp2 = fix(rand(1,1)*numGroups)+1;
    end
    
    % Switch group
    newGroups = curGroups;
    newGroups(ind) = grp2;
    
    % Evaluate likelihood
    newlik = 0;
    for i=1:(numEvents-1)
        for j=(i+1):numEvents
            if(newGroups(i)<newGroups(j))
                newlik = newlik + logp(i,j)+ logmp(j,i);
            elseif(newGroups(i)==newGroups(j))
                newlik = newlik + logmp(i,j) + logmp(j,i);
            else
                newlik = newlik + logp(j,i) + logmp(i,j);
            end
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
        curGroups = newGroups;
    end

    if(it>burnin && mod(it-burnin, interval)==0)
        samples(sno,:) = curGroups;
        liks(sno) = curlik;
        sno = sno+1;
        display(sprintf('It: %i. Lik: %f', it, curlik));
    end
end



%% Look at the likelihood distribution
figure; hist(liks, 100);
%xlim([max(liks)-50, max(liks)]);




%% Extract just the unique samples in ascending order of likelihood.
[uliks ulikinds n] = unique(liks);
samples(ulikinds,:)

% Use the maximum likelihood grouping
groupLabels = samples(ulikinds(end), :);



%% Compute the likelihood of this number of clusters
% Marginalize of the cluster assignment
pgroupnum = mean(exp(liks))
pngrp(numGroups) = pgroupnum;
t = [mean(liks) max(liks)]
likngrp(numGroups, :) = t;



%% Save the results
%save PatientsBothHemiSeqGroupsOutput.mat samples burnin interval numSamples  p groupLabels atrophy
%save PatientsLeftHemiSeqGroupsOutput.mat samples burnin interval numSamples p groupLabels atrophy
%save PatientsRightHemiSeqGroupsOutput.mat samples burnin interval numSamples  p groupLabels atrophy
%save ControlsBothHemiSeqGroupsOutput.mat samples burnin interval numSamples  p groupLabels atrophy
%save ControlsLeftHemiSeqGroupsOutput.mat samples burnin interval numSamples  p groupLabels atrophy
%save ControlsRightHemiSeqGroupsOutput.mat samples burnin interval numSamples  p groupLabels atrophy

save PatientsBothHemiSeqGroups08_P2_Output.mat samples burnin interval numSamples groupLabels p PatAtrophy ConAtrophy patjdiff



%% Visualize the maximum likelihood grouping.

figure; 

% Compute the within group likelihoods and normalize by number of group members.
for i=1:numGroups
    groupinds = find(groupLabels==i);
    ngrouplik(i) = 0;
    grpsize = length(groupinds);
    for j=1:(grpsize-1)
        for k=(j+1):grpsize
            ngrouplik(i) = ngrouplik(i)+log(p(groupinds(j), groupinds(k)))+log(p(groupinds(k), groupinds(j)));
        end
    end
    ngrouplik(i) = ngrouplik(i)/(grpsize*(grpsize+1)/2);
end

% Compute normalized predecession likelihoods for consecutive groups.
npredlik(1) = 0;
for i=1:(numGroups-1)
    groupinds = find(groupLabels==i);
    ngroupinds = find(groupLabels==(i+1));
    npredlik(i+1) = 0;
    for j=1:(length(groupinds))
        for k=1:length(ngroupinds)
            npredlik(i+1) = npredlik(i+1)+log(p(groupinds(j), ngroupinds(k)));
        end
    end
    npredlik(i+1) = npredlik(i+1)/(length(groupinds)*length(ngroupinds));
end



% Within each group, the vertical positions vary, but horizontal positions
% are the same.
vertPos = zeros(length(groupLabels), 1);
horPos = zeros(length(groupLabels), 1);
colmax = 9;
for i=1:length(groupLabels)
    grpsize = sum(groupLabels(1:(i-1))==groupLabels(i));
    vertPos(i) = mod(grpsize,colmax);
%    horPos(i) = sum(exp(npredlik(1:groupLabels(i))));
    horPos(i) = groupLabels(i) + fix(grpsize/colmax)*0.3;
end

for i=1:35
    labels{i} = sprintf('%02iL', i);
    labels{i+35} = sprintf('%02iR', i);
end

scatter(horPos, vertPos, 800*ones(length(groupLabels),1), patjdiff, 'SizeData', 800, 'LineWidth', 3);
colormap winter
hold on;
for i=1:numEvents;
%    text(horPos(i)-0.15, vertPos(i), labels{i}(1:2));
    text(horPos(i)-0.1, vertPos(i), labels{i});
end
xlim([0 numGroups+1]);
ylim([-1 max(vertPos)+1]);

%title(sprintf('%i stages of atrophy (both hemispheres)', numGroups));
%title(sprintf('%i stages of atrophy (left hemisphere)', numGroups));
%title(sprintf('%i stages of atrophy (right hemisphere)', numGroups));
axis off;

%% Write out the figure
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 14 6]);
%print -dpng 'PatientsLeftHemiSeqGroups.png'
%print -dpng 'PatientsRightHemiSeqGroups.png'
print -dpng 'PatientsBothHemiSeqGroupsP2.png'
%print -dpng 'ControlsLeftHemiSeqGroups.png'
%print -dpng 'ControlsRightHemiSeqGroups.png'
%print -dpng 'ControlsBothHemiSeqGroups.png'

