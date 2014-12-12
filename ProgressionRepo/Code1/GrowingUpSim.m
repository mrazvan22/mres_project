%% Simulation of growing-up data.

%% Definition of probability distribution on the date of various events
% 11 +/- 1.5 year; min 6 years
pubonset = @(t)(max([(randn(1,t)*1.5 + 11)*365.25; ones(1,t)*6*365.25]));

% 6y +/- 6 m
adultteeth = @(t)((randn(1,t)*6/12 + 6)*365.25);

% 13 +/- 3 month; min 4 months.
firstword = @(t)(max([(randn(1,t)*3 + 13)*365.25/12; ones(1,t)*4*365.25/12]));

% 13 +/- 2 year; min 5 years
firstkiss = @(t)(max([(randn(1,t)*2 + 13)*365.25; ones(1,t)*5*365.25]));

% Exactly at 13y
fourteenthbday = @(t)(ones(1,t)*13*365.25);

% 7 +/- 3 month; min 1 month
firsttooth = @(t)(max([(randn(1,t)*3 + 7)*365.25/12; ones(1,t)*365.25/12]));

% uniform between 4 and 5 years
startschool = @(t)((rand(1,t) + 4)*365.25);


%% Create a set of event dates for a number of individuals and plot distributions.

samples = 1000;
dates = [pubonset(samples); adultteeth(samples); firstword(samples); firstkiss(samples); fourteenthbday(samples); firsttooth(samples); startschool(samples)];

% Plot histograms of each
figure;
hold on;
[f1 l1] = hist(dates(1,:),100);
bar(l1, f1, 'b');
[f2 l2] = hist(dates(2,:),20);
bar(l2, f2, 'r');
[f3 l3] = hist(dates(3,:),20);
bar(l3, f3, 'g');
[f4 l4] = hist(dates(4,:),100);
bar(l4, f4, 'm');
[f5 l5] = hist(dates(5,:),[13*365.25-50 13*365.25 13*365.25+50]);
bar(l5, f5, 'c');
[f6 l6] = hist(dates(6,:),20);
bar(l6, f6, 'y');
[f7 l7] = hist(dates(7,:),20);
bar(l7, f7, 'k');
ylim([0 200]);
legend('Puberty', 'Adult Teeth', 'First Word', 'First Kiss', '13th Bday', 'First Tooth', 'Start School');
xlabel('Age/days');
ylabel('Bin count');

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 12 9]);
print('-dpng', 'GrowingUpDists.png')


%% Construct precedence matrix.
[numEvents numSamples] = size(dates);
p1 = ones(numEvents, numEvents);
p2 = ones(numEvents, numEvents);
for i=1:(numEvents-1)
    for j=(i+1):numEvents
        ipj = find(dates(i,:)<dates(j,:));
        p1(i,j) = length(ipj)/numSamples;
        p1(j,i) = 1-p1(i,j);
        p2(i,j) = sum(dates(j,ipj) - dates(i,ipj))/sum(abs(dates(i,:) - dates(j,:)));
        p2(j,i) = 1-p2(i,j);
    end
end

figure;
subplot(2,2,1);
imagesc(p1);
subplot(2,2,2);
imagesc(p2);

figure;
% Add matrices for increasingly small data sets
ind = 0;
for s = [1000 500 100 50 10 5];
    subdates = dates(:,1:s);
    p1 = ones(numEvents, numEvents);
    p2 = ones(numEvents, numEvents);
    for i=1:(numEvents-1)
        for j=(i+1):numEvents
            ipj = find(subdates(i,:)<subdates(j,:));
            p1(i,j) = length(ipj)/s;
            p1(j,i) = 1-p1(i,j);
            p2(i,j) = sum(subdates(j,ipj) - subdates(i,ipj))/sum(abs(subdates(i,:) - subdates(j,:)));
            p2(j,i) = 1-p2(i,j);
        end
    end
    ind = ind + 1;
    subplot(6,2,ind);
    imagesc(p1);
    if(ind == 1)
        title('Summed');
    end
    %axis off;
    ylabel(sprintf('N = %i', s));
    ind = ind + 1;
    subplot(6,2,ind);
    imagesc(p2);
    if(ind == 2)
        title('Weighted');
    end
    %axis off;
end

% p2 is more robust to small numbers of samples, but may be more vulnerable
% to outliers or errors in the input data.

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6 9]);
print('-dpng', 'GU_DirRetroPs.png')



%% Now construct matrices from snapshot data.

% Create an age for each subject
[numEvents numSamples] = size(dates);
ages = rand(1,numSamples)*max(max(dates));
knowndates = dates;
knowndates(find(repmat(ages, [numEvents 1])<dates)) = 0;

figure;
hold on;
ind = 0;
for s = [1000 500 100 50 10 5];
    subdates = knowndates(:,1:s);
    p1 = ones(numEvents, numEvents);
    p2 = ones(numEvents, numEvents);
    p3 = ones(numEvents, numEvents);
    for i=1:(numEvents-1)
        for j=(i+1):numEvents
            ixorj = find(xor(subdates(i,:) == 0, subdates(j,:) == 0));
            inotj = find(subdates(i,:) > 0 & subdates(j,:) == 0);
            jnoti = find(subdates(j,:) > 0 & subdates(i,:) == 0);
            p1(i,j) = length(inotj)/length(ixorj);
            p1(j,i) = 1-p1(i,j);
            p2(i,j) = sum(ages(inotj) - subdates(i,inotj))/(sum(ages(inotj) - subdates(i,inotj)) + sum(ages(jnoti) - subdates(j,jnoti)));
            p2(j,i) = 1-p2(i,j);
            sdi = subdates(i,:);
            sdj = subdates(j,:);
            sdi(find(sdi==0)) = ages(find(sdi==0));
            sdj(find(sdj==0)) = ages(find(sdj==0));
            ipj = find(sdi<sdj);
            p3(i,j) = sum(sdj(ipj) - sdi(ipj))/sum(abs(sdi - sdj));
            p3(j,i) = 1-p3(i,j);
        end
    end
    ind = ind + 1;
    subplot(6,3,ind);
    imagesc(p1);
    if(ind == 1)
        title('Summed');
    end
    ylabel(sprintf('N = %i', s));
    ind = ind + 1;
    subplot(6,3,ind);
    imagesc(p2);
    if(ind == 2)
        title('Weighted');
    end
    ind = ind + 1;
    subplot(6,3,ind);
    imagesc(p3);
    if(ind == 3)
        title('Inclusive');
    end
end

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 9 9]);
print('-dpng', 'GU_DirSnapshotPs.png')


