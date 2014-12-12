%% Simulation of normal and AD atrophy for full simulation of ordering process.

% Relative volume of whole brain.
meanbrainvol = 1;
stdbrainvol = 0.2;

% Mean relative volumes of each region
meanhipvol = 0.1;
stdhipvol = meanhipvol/10;
meanprecvol = 0.05;
stdprecvol = meanprecvol/10;
meantlcvol = 0.15;
stdtlcvol = meantlcvol/10;
meanpcgvol = 0.06;
stdpcgvol = meanpcgvol/10;
meanaccvol = 0.02;
stdaccvol = meanaccvol/10;

NumPatients = 1000;
NumControls = 1000;

% Create regional volumes before onset for each patient and control
PatBrainVols = randn(1,NumPatients)*stdbrainvol + meanbrainvol;
ConBrainVols = randn(1,NumControls)*stdbrainvol + meanbrainvol;

PatInitRegRelVols = randn(NumPatients, 5).*repmat([stdhipvol, stdprecvol, stdtlcvol stdpcgvol stdaccvol], [NumPatients 1])...
    +repmat([meanhipvol, meanprecvol, meantlcvol meanpcgvol meanaccvol], [NumPatients 1]);
PatInitRegVols = repmat(PatBrainVols, [5,1])'.*PatInitRegRelVols;
ConInitRegRelVols = randn(NumPatients, 5).*repmat([stdhipvol, stdprecvol, stdtlcvol stdpcgvol stdaccvol], [NumControls 1])...
    +repmat([meanhipvol, meanprecvol, meantlcvol meanpcgvol meanaccvol], [NumControls 1]);
ConInitRegVols = repmat(ConBrainVols, [5,1])'.*ConInitRegRelVols;

% Disease duration in years.
meanduration = 15;
stdduration = 3;

PatDurations = randn(1,NumPatients)*stdduration + meanduration;
maxdur = max(PatDurations);
% Sample volume monthly.
SampleTimes = 0:(1/12):maxdur;

% Mean and std of onset times in different regions.  These are fractions of disease
% duration.
meanhiponset = 0.1;
stdhiponset = 0.02;
meanpreconset = 0.2;
stdpreconset = 0.02;
meantlconset = 0.4;
stdtlconset = 0.02;
meanpcgonset = 0.4;
stdpcgonset = 0.02;
meanacconset = 0.7;
stdacconset = 0.02;

% Onset of atrophy in each region varies but linked to overall disease
% duration.
HipOnsets = (randn(1,NumPatients)*stdhiponset + meanhiponset).*PatDurations;
HipOnsets(find(HipOnsets<=0)) = min(HipOnsets(find(HipOnsets>0)));
PrecOnsets = (randn(1,NumPatients)*stdpreconset + meanpreconset).*PatDurations;
PrecOnsets(find(PrecOnsets<=0)) = min(PrecOnsets(find(PrecOnsets>0)));
TLC_Onsets = (randn(1,NumPatients)*stdtlconset + meantlconset).*PatDurations;
TLC_Onsets(find(TLC_Onsets<=0)) = min(TLC_Onsets(find(TLC_Onsets>0)));
PCG_Onsets = (randn(1,NumPatients)*stdpcgonset + meanpcgonset).*PatDurations;
PCG_Onsets(find(PCG_Onsets<=0)) = min(PCG_Onsets(find(PCG_Onsets>0)));
AccOnsets = (randn(1,NumPatients)*stdacconset + meanacconset).*PatDurations;
AccOnsets(find(AccOnsets<=0)) = min(AccOnsets(find(AccOnsets>0)));

% Duration of atrophy in each region varies but linked to overall disease
% duration.
HipDuration = randn(1,NumPatients)*stdduration + PatDurations;
HipDuration(find(HipDuration<=0)) = min(HipDuration(find(HipDuration>0)));
PrecDuration = randn(1,NumPatients)*stdduration + PatDurations;
PrecDuration(find(PrecDuration<=0)) = min(PrecDuration(find(PrecDuration>0)));
TLC_Duration = randn(1,NumPatients)*stdduration + PatDurations;
TLC_Duration(find(TLC_Duration<=0)) = min(TLC_Duration(find(TLC_Duration>0)));
PCG_Duration = randn(1,NumPatients)*stdduration + PatDurations;
PCG_Duration(find(PCG_Duration<=0)) = min(PCG_Duration(find(PCG_Duration>0)));
AccDuration = randn(1,NumPatients)*stdduration + PatDurations;
AccDuration(find(AccDuration<=0)) = min(AccDuration(find(AccDuration>0)));

% Minimum regional volume as a fraction of the initial volume.
minvolprop = 0.5;

% Create volumes relative to initial volume using a cumulative normal
% distribution model.  At onset time, the argument to the normcdf is -3 and
% at the end of the atrophy duration it is +3.
relhipvols = (1- minvolprop*normcdf((6*repmat(SampleTimes, [NumPatients 1]) - 6*repmat(HipOnsets', [1 length(SampleTimes)]) - 3*repmat(HipDuration', [1 length(SampleTimes)]))./repmat(HipDuration', [1 length(SampleTimes)])));
relhipvols(find(relhipvols<=0)) = min(relhipvols(find(relhipvols>0)));
relprecvols = (1- minvolprop*normcdf((6*repmat(SampleTimes, [NumPatients 1]) - 6*repmat(PrecOnsets', [1 length(SampleTimes)]) - 3*repmat(PrecDuration', [1 length(SampleTimes)]))./repmat(PrecDuration', [1 length(SampleTimes)])));
relprecvols(find(relprecvols<=0)) = min(relprecvols(find(relprecvols>0)));
reltlcvols = (1- minvolprop*normcdf((6*repmat(SampleTimes, [NumPatients 1]) - 6*repmat(TLC_Onsets', [1 length(SampleTimes)]) - 3*repmat(TLC_Duration', [1 length(SampleTimes)]))./repmat(TLC_Duration', [1 length(SampleTimes)])));
reltlcvols(find(reltlcvols<=0)) = min(reltlcvols(find(reltlcvols>0)));
relpcgvols = (1- minvolprop*normcdf((6*repmat(SampleTimes, [NumPatients 1]) - 6*repmat(PCG_Onsets', [1 length(SampleTimes)]) - 3*repmat(PCG_Duration', [1 length(SampleTimes)]))./repmat(PCG_Duration', [1 length(SampleTimes)])));
relpcgvols(find(relpcgvols<=0)) = min(relpcgvols(find(relpcgvols>0)));
relaccvols = (1- minvolprop*normcdf((6*repmat(SampleTimes, [NumPatients 1]) - 6*repmat(AccOnsets', [1 length(SampleTimes)]) - 3*repmat(AccDuration', [1 length(SampleTimes)]))./repmat(AccDuration', [1 length(SampleTimes)])));
relaccvols(find(relaccvols<=0)) = min(relaccvols(find(relaccvols>0)));

HipVols = relhipvols.*repmat(PatInitRegVols(:,1), [1 length(SampleTimes)]);
PrecVols = relprecvols.*repmat(PatInitRegVols(:,2), [1 length(SampleTimes)]);
TLC_Vols = reltlcvols.*repmat(PatInitRegVols(:,3), [1 length(SampleTimes)]);
PCG_Vols = relpcgvols.*repmat(PatInitRegVols(:,4), [1 length(SampleTimes)]);
AccVols = relaccvols.*repmat(PatInitRegVols(:,5), [1 length(SampleTimes)]);

% Check out an example.
figure; plot(HipVols(2,:)); hold on; plot(PrecVols(2,:)); plot(TLC_Vols(2,:)); plot(PCG_Vols(2,:)); plot(AccVols(2,:));


%% Ground truth P matrix from the actual onset times.
numEvents = 5;
pgt1 = ones(numEvents, numEvents);
pgt2 = ones(numEvents, numEvents);
onsets = [HipOnsets; PrecOnsets; TLC_Onsets; PCG_Onsets; AccOnsets];
for i=1:(numEvents-1)
    for j=(i+1):numEvents
        ipj = find(onsets(i,:)<onsets(j,:));
        pgt1(i,j) = length(ipj)/NumPatients;
        pgt1(j,i) = 1-pgt1(i,j);
        pgt2(i,j) = sum(onsets(j,ipj) - onsets(i,ipj))/sum(abs(onsets(i,:) - onsets(j,:)));
        pgt2(j,i) = 1-pgt2(i,j);
    end
end


%% Estimate P from all time points using known initial volumes.

% Need to assume a certain variance on volume measurements to obtain a
% probability that the volume has been reduced.
volmeasstd = 0.02;
% Now we can compute the probability that the measurement at time t came
% from the same true volume as time 0.
AllVols(1,:,:) = HipVols;
AllVols(2,:,:) = PrecVols;
AllVols(3,:,:) = TLC_Vols;
AllVols(4,:,:) = PCG_Vols;
AllVols(5,:,:) = AccVols;

probmeas = zeros(numEvents, NumPatients, length(SampleTimes));
for i=1:numEvents
    for j=1:NumPatients
        for k=1:length(SampleTimes)
            Z = (AllVols(i,j,1) - AllVols(i,j,k))/volmeasstd;
            probmeas(i,j,k) = 2 - 2*cdf('norm', Z, 0, 1);
        end
    end
end

probE = 1-probmeas;

%figure; imagesc(squeeze(probE(1,:,:)))
%figure; imagesc(squeeze(probE(2,:,:)))
%figure; imagesc(squeeze(probE(3,:,:)))
%figure; imagesc(squeeze(probE(4,:,:)))
%figure; imagesc(squeeze(probE(5,:,:)))

%figure; imshow(squeeze(probE(1,:,:).*(1-probE(2,:,:))),[0 1])
%figure; imshow(squeeze(probE(2,:,:).*(1-probE(1,:,:))),[0 1])

p1 = ones(numEvents, numEvents);
for i=1:(numEvents-1)
    for j=(i+1):numEvents
        p1(i,j) = sum(sum((squeeze(probE(i,:,:)).*(1-squeeze(probE(j,:,:)))).^2))/sum(sum(squeeze(probE(i,:,:)).*(1-squeeze(probE(j,:,:)))));
        p1(j,i) = 1-p1(i,j);
    end
end

% Compute a separate p for each subject
for i=1:NumPatients
    for j=1:(numEvents-1)
        for k=(j+1):numEvents
            p2ind(i,j,k) = max(probE(j,i,:).*(1-probE(k,i,:)));
            p2ind(i,k,j) = max(probE(k,i,:).*(1-probE(j,i,:)));
        end
    end
end
p2a = squeeze(mean(p2ind));
p2b = squeeze(sum(p2ind.^2))./max(squeeze(sum(p2ind)),1E-10);

% p1 doesn't really work.  p2 seems OK; noticably not antisymmetric, but that is
% not necessarily a problem.


%% Now from one random time point from each patient with initial volumes known

SnapShotIndices = fix(rand(1,NumPatients)*length(SampleTimes))+1;
SnapShotVols = zeros(numEvents, NumPatients);
for i=1:NumPatients
    for j=1:numEvents
        SnapShotVols(j,i) = AllVols(j,i,SnapShotIndices(i));
    end
end

probmeas = zeros(numEvents, NumPatients);
for i=1:numEvents
    for j=1:NumPatients
        Z = (AllVols(i,j,1) - SnapShotVols(i,j))/volmeasstd;
        probmeas(i,j) = 2 - 2*cdf('norm', Z, 0, 1);
    end
end
probE = 1-probmeas;

p3 = ((probE.^2)*((1-probE').^2))./(probE*(1-probE'));


%% Now without knowledge of the initial volume, but using statistics from controls.

NormRegVols = SnapShotVols./repmat(PatBrainVols, [numEvents 1]);

% Now base the probability of atrophy having started on whether the normalized volume
% of the region is outside the normal range.
NormControlRegVols = ConInitRegVols./repmat(ConBrainVols, [numEvents 1])';
meanregnormvols = mean(NormControlRegVols);
stdregnormvols = std(NormControlRegVols);

probmeas = zeros(numEvents, NumPatients);
for i=1:numEvents
    for j=1:NumPatients
        Z = (meanregnormvols(i) - NormRegVols(i,j))/stdregnormvols(i);
        if(Z<0)
            Z=0;
        end
        probmeas(i,j) = 2 - 2*cdf('norm', Z, 0, 1);
    end
end
probE = 1-probmeas;

p4 = ((probE.^2)*((1-probE').^2))./(probE*(1-probE'));





