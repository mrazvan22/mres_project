%% Loading data
clear all
close all

file_data = spm_select(1, 'mat', 'Give original data');
file_order = spm_select(1, 'mat', 'Give orderings');
load(file_data);
load(file_order);

nr_pat = size(data_struct.data_patient{1}.data_nomixt, 2);
nr_roi = length(data_struct.name_roi);

%% Postprocessing

avPos = mean(order_events_thinned, 2);
stdPos = std(order_events_thinned, [], 2);
[ord inds_av] = sort(avPos);
[ord2, inds2] = sort(inds_av);

%% Producing goose plots

% Horizontal positions of each marker are just their average position in
% the order.
% Vertical positions are displaced to avoid overlap of markers and
% emphasize the clustering.  If a marker does not overlap with its
% predecessor, place it at vertical position of zero.  Otherwise move it up
% far enough not to overlap.  For sequences of overlapping ones, we move
% them alternately up and down.

% For 34 regions
%xRad = 0.5;
%yRad = 0.04;
% For 68 regions
rx = 0.4;
ry = 0.4;
vertCo = zeros(size(avPos));
labels{1} = sprintf('%02i', 1);
labels{2} = sprintf('%02iR', 1);
for i=1:35
    labels{i+2} = sprintf('%02i', i+1);
    labels{i+36} = sprintf('%02iR', i+1);
end
labels_diagnosis = 'ABC';

figure(2), clf
t = 0:pi/120:2*pi;
px = zeros(length(t), 1);
py = zeros(length(t), 1);
colourRGB = [30 144 255]/255;
colourRGBDiag = [255 127 80]/255;
vertCo = 1;
for roi = 1:nr_roi,
    
    px = stdPos(inds_av(roi))*cos(t) + avPos(inds_av(roi));
    py = 1*sin(t) + vertCo;
    pp = patch(px, py, colourRGB, 'EdgeColor', [0 0 1], 'LineWidth', 2);
    text(avPos(inds_av(roi)), vertCo, sprintf('%02i', inds_av(roi)), ...
        'HorizontalAlignment', 'center', 'Color', [255 255 255]/255, 'FontSize', 12, 'FontWeight', 'demi');
    vertCo = vertCo + 2;

end
hold on
axis equal, axis off

vertCo = zeros(size(avPos));
for i=2:length(avPos),
    
    if find([6:5:numEvents] == i),
        
        vertCo(inds_av(i)) = 0;
        
    else
        
        vertCo(inds_av(i)) = (vertCo(inds_av(i-1)) + 3*ry);
        
    end
    
end


figure(3), clf
stdPosFull(stdPosFull == 0) = 0.5;
t = 0:pi/120:2*pi;
px = zeros(length(t), 1);
py = zeros(length(t), 1);
colourRGB = [100 149 237]/255;
colourRGBDiag = [255 127 80]/255;
plot(avPos(inds_av), vertCo(inds_av), '--k')
for i = 1:numEvents,
    
    px = stdPosFull(i)*cos(t) + avPos(i);
    py = ry*sin(t) + vertCo(i);
    
    if i <= numRegions,
        
        pp(i) = patch(px, py, colourRGB, 'EdgeColor', [0 0 1], 'LineWidth', 2);
        text(avPos(i), vertCo(i), sprintf('%02i', i), ...
            'HorizontalAlignment', 'center', 'Color', [255 255 255]/255, 'FontSize', 12, 'FontWeight', 'demi');
        
    else
        
        pp(i) = patch(px, py, colourRGBDiag, 'EdgeColor', [1 0 0], 'LineWidth', 2);
        text(avPos(i), vertCo(i), labels_diagnosis(i-numRegions), ...
            'HorizontalAlignment', 'center', 'Color', [255 255 255]/255, 'FontSize', 12, 'FontWeight', 'demi');
        
    end
    
end
hold on
axis equal, axis off

% Write out the figure
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 14 6]);

%% Making an occurrence histogram...
hist2_mat = zeros(nr_roi, nr_roi);
for roi1 = 1:nr_roi,
    
    for roi2 = 1:nr_roi,
        
        hist2_mat(roi1, roi2) = sum(order_events_thinned(inds_av(roi1), :) == roi2);
        
    end
    
end
figure(3), clf
imagesc(log(hist2_mat))


%% Studying individual order of patients

patient_id = data_struct.patient_id;
patient_tp = data_struct.patient_tp;
model_atrophy = zeros(numEvents, numEvents);
for i = 1:numEvents,
    
    model_atrophy(i, round(avPos(i)):end) = 1;
    
end

for tp = 1:numMeas,
    
    atrophy_local = atrophy_sig(:, tp);
    atrophy_local = repmat(atrophy_local, 1, numEvents);
    overlap_mat_ones = sum((atrophy_local == 1) &(model_atrophy == 1));
    overlap_mat_zeros = sum((atrophy_local == 0) &(model_atrophy == 0));
    overlap_mat = (overlap_mat_ones + overlap_mat_zeros)/numEvents;
    I = find(overlap_mat == max(overlap_mat));
    if length(I) == 1,
        
        class_tp(tp) = I;
        overlap_tp(tp) = max(overlap_mat);
        
    else
        
        I2 = find(overlap_mat_ones(I) == max(overlap_mat_ones(I)));
        class_tp(tp) = I(I2(1));
        overlap_tp(tp) = overlap_mat(I(I2(1)));
        
    end
    
end


for id = unique(patient_id),
    
    I_pt = find(patient_id == id);
    patient_tp(I_pt)
    class_tp(I_pt),
    pause
    
end

rx = 1;
ry = 1;

vertCo = zeros(size(avPos));
for i=2:length(avPos),
    
    if find([6:5:38] == i),
        
        vertCo(inds_av(i)) = 0;
        
    elseif(vertCo(inds_av(i-1))>=0)
        
        vertCo(inds_av(i)) = -(vertCo(inds_av(i-1)) + 2*ry);
        
    else
        
        vertCo(inds_av(i)) = -vertCo(inds_av(i-1));
        
    end
    
end

figure(4), clf
t = 0:pi/120:2*pi;
px = zeros(length(t), 1);
py = zeros(length(t), 1);
colourRGB = [30 144 255]/255;
colourRGBDiag = [255 127 80]/255;
for i = 1:numEvents,
    
    px = rx*cos(t) + inds2(i);
    py = ry*sin(t) + vertCo(i);
    
    if i <= numRegions,
        
        pp(i) = patch(px, py, colourRGB, 'EdgeColor', [0 0 1], 'LineWidth', 2);
        text(inds2(i), vertCo(i), sprintf('%02i', i), ...
            'HorizontalAlignment', 'center', 'Color', [255 255 255]/255, 'FontSize', 12, 'FontWeight', 'demi');
        
    else
        
        pp(i) = patch(px, py, colourRGBDiag, 'EdgeColor', [1 0 0], 'LineWidth', 2);
        text(inds2(i), vertCo(i), labels_diagnosis(i-numRegions), ...
            'HorizontalAlignment', 'center', 'Color', [255 255 255]/255, ...
            'FontSize', 10, 'FontWeight', 'demi', 'FontName', 'AvantGarde');
        
    end
    
end
hold on

maxVertCo = max(vertCo);
vertCo = min(vertCo);
for id = unique(patient_id),
    
    colourRGB = rand(1, 3);
    vertCo = vertCo - 4*ry;
    I_pt = find(patient_id == id);
    tp = patient_tp(I_pt);
    class_pt = class_tp(I_pt);
    vertCo_local = repmat(vertCo, 1, length(tp));
    class_unique = unique(class_pt);
    for c = 1:length(class_unique),
        
       I = find(class_pt == class_unique(c));
       if length(I) > 1,
                    
           vertCo_adjust = ry*linspace(-1, 1, length(I));
           vertCo_local(I) = vertCo_local(I) + vertCo_adjust;
           
       end
       
    end
                     
    for i = 1:length(I_pt),

        px = rx*cos(t) + class_pt(i);
        py = ry*sin(t) + vertCo_local(i);

        patch(px, py, colourRGB, 'EdgeColor', colourRGB, 'LineWidth', 2);
        text(class_pt(i), vertCo_local(i), sprintf('t%2d', tp(i)), ...
            'HorizontalAlignment', 'center', 'Color', [255 255 255]/255, ...
            'FontSize', 10, 'FontWeight', 'demi', 'FontName', 'AvantGarde');
               
    end
    
end
axis equal, axis off        

%% Making color maps

map = colormap('hot');
A = zeros(numEvents, 4);
A(:, 1) = 1:numEvents;
for i = 1:3,
    
    A(:, i+1) = resample(map(:, i), numEvents, size(map, 1));
    A(A(:, i+1) < 0, i+1) = 0;
    A(A(:, i+1) > 1, i+1) = 1;
    A(:, i+1) = round(255*A(:, i+1));
    
end

for i = 1:numEvents,
    
    fprintf('%d %d %d %d 1 1 1 "Label %d"\n', ...
        A(i, 1), A(i, 2), A(i, 3), A(i, 4), A(i, 1))
    
end
