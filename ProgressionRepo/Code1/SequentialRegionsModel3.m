clear all
close all

file_data = spm_select(1, 'mat');
load(file_data);
%% Create the precedence probability matrix.
atrophy_patient = data_struct.data_jacc_patient_select';
atrophy_control = data_struct.data_jacc_control_select';
I_lh = [2 3:36];
I_rh = [1 37:size(atrophy_patient, 2)];
atrophy_patient_lhrh(:, :, 1) = atrophy_patient(:, I_lh);
atrophy_patient_lhrh(:, :, 2) = atrophy_patient(:, I_rh);
atrophy_patient_mean = mean(atrophy_patient_lhrh, 3);
atrophy_control_lhrh(:, :, 1) = atrophy_control(:, I_lh);
atrophy_control_lhrh(:, :, 2) = atrophy_control(:, I_rh);
atrophy_control_mean = mean(atrophy_control_lhrh, 3);

I_label_select = data_struct.I_label_select;

[numMeas numEvents] = size(atrophy_patient_mean);
p1 = PrecedenceMatrix1(atrophy_patient_mean);
p2 = PrecedenceMatrix2(atrophy_patient_mean, atrophy_control_mean);

p = p2;
p(p == 0) = eps;
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

burnin = 1e5;
interval = 500;
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

order_events = zeros(numSamples, numEvents);
for it = 1:numSamples,
    
    [d, order_events(it, :)] = sort(samples(it, :));
    
end
    
hist2_mat = zeros(numEvents);
[ord inds_av] = sort(mean(order_events(:, :)));
for e1 = 1:numEvents,
    
    for e2 = 1:numEvents,
        
        hist2_mat(e1, e2) = sum(order_events(:, inds_av(e1)) == e2);
        
    end
    
end
[ord, sort_inds_av] = sort(inds_av);

name_roi{1} = data_struct.name_roi{1};
for e = 2:(numEvents+1),
    
    name_roi{e} = data_struct.name_roi{e+1};
    
end

max_length_str = 0;
for e = 1:numEvents,
    
    max_length_str = max([max_length_str length(name_roi{e})]);
    
end
mat_labels = zeros(numEvents, max_length_str);
for e = 1:numEvents,
    
    length_str = length(name_roi{inds_av(e)});
    mat_labels(e, 1:length_str) = name_roi{inds_av(e)};
    
end
figure(1), clf
imagesc(hist2_mat(numEvents:-1:1, numEvents:-1:1))
map_inverted = repmat(linspace(1, 0, 64)', 1, 3);
colormap(map_inverted)
axis square
set(gca, ...
    'YTick', [1:numEvents], 'YTickLabel', char(mat_labels(numEvents:-1:1, :)), ...
    'YGrid', 'on')
hXLabel = xlabel('Model place');
hYLabel = ylabel('Region');
set([hXLabel hYLabel], ...
    'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'demi');


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



%% Compute the average position of each event weighted by path probabilities
sumord=zeros(1,numEvents);
sumprobs = 0;
for i=1:length(ulikinds)
    [a b] = sort(samples(ulikinds(i),:));
    prob = exp(uliks(i));
    sumord = sumord+b*prob;
    sumprobs = sumprobs + prob;
end
avPos = round(sumord/sumprobs)
order_reorder = [35:-1:1];
for e = 1:numEvents,
    
    avPos_reorder(avPos == e) = order_reorder(e);
    
end
    

file_segm = 'img_segm.nii';
V_segm = spm_vol(file_segm);
img_segm = spm_read_vols(V_segm);
img_order = zeros(V_segm.dim);
for e = 1:numEvents,
    
    img_order(I_label_select{I_lh(e)}) = avPos_reorder(e);
    
end
V_order = V_segm;
V_order.fname = 'img_order_data_AD_precmat2.nii';
spm_create_vol(V_order);
spm_write_vol(V_order, img_order);
%% Save the results
save PatientsBothHemiSeqRegionsOutput.mat samples burnin interval numSamples avPos p atrophy
save ControlsBothHemiSeqRegionsOutput.mat samples burnin interval numSamples avPos p atrophy


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
% For 34 regions
%xRad = 0.5;
%yRad = 0.04;
% For 68 regions
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

labels{1} = sprintf('%02iL', 1);
labels{2} = sprintf('%02iR', 1);
for i=1:34
    labels{i+2} = sprintf('%02iL', i+1);
    labels{i+36} = sprintf('%02iR', i+1);
end
% One hemisphere
scatter(avPos, vertCo, 'SizeData', 500, 'LineWidth', 3);
% Both hemispheres
scatter(avPos, vertCo, 'SizeData', 1000, 'LineWidth', 3);
hold on;
for i=1:length(avPos);
    % This one is for both hemispheres (68 or 70 regions)
    text(avPos(i)-0.7, vertCo(i), labels{i});
    % This one works for the 34 or 35 regions (single hemisphere)
    %text(avPos(i)-0.2, vertCo(i), sprintf('%02i', i));
end
% One hemi
%ylim([min(vertCo),-min(vertCo)]*5/2);
% Both hemis
ylim([min(vertCo),-min(vertCo)]);

axis off;

%% Write out the figure
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 14 6]);
%print -dpng 'fAD_34_LH_Ordering.png'
%print -dpng 'fAD_34_RH_Ordering.png'
%print -dpng 'fAD_68_BothH_Ordering.png'
print -dpng 'PatientsBothHemiSeqRegionsOrdering.png'
print -dpng 'ControlsBothHemiSeqRegionsOrdering.png'



%% Alternative very simple ordering from average atrophy in each region
avPos = mean(atrophy);
avPos = (avPos-min(avPos)/(max(avPos)-min(avPos));
avPos = (avPos-min(avPos))/(max(avPos)-min(avPos));
avPos = (avPos-min(avPos))/(max(avPos)-min(avPos))*68
% This produces a similar ordering to the above.


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
scatter(AvPosPatient([1 3:36]), AvPosPatient([2 37:end]), 'SizeData', 500, 'LineWidth', 3);
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


%% Compute yes or no atrophy matrix based on comparison of patients with control distribution.
load data_jacob_control.mat
ConAtrophy = data_jacc_control([(end-69):end], :)';
mat = mean(ConAtrophy);
sat = std(ConAtrophy);

load data_jacob_patient.mat
PatAtrophy = data_jacc_patient([(end-69):end], :)';
[numMeas numEvents] = size(PatAtrophy);

pat = zeros(size(PatAtrophy));
for i=1:numMeas
    for j=1:numEvents
        pat(i,j) = exp(-(PatAtrophy(i,j) - mat(i))^2/(2*sat(i)^2))/(sqrt(2*pi*sat(i)^2));
    end
end

