%plotKsDensityPerRegion

%load('Mat Data/MAPER_MASKED_Data.mat');
%load('/Users/jaghajan/Documents/Experiments/Mat Data/MAPER_Data.mat');
load('/Users/jaghajan/Documents/Experiments/Mat Data/MAPER_ICV_NORM_DATA.mat')

mci_indx=find(data_struct.status_mat==2);
ad_indx=find(data_struct.status_mat==3);
nr_roi=83;

for roi=1:nr_roi
    [f_controls,x_controls]=ksdensity(squeeze(data_struct.data_control.avg_vol(roi,:)));
    [f_mci,x_mci]=ksdensity(squeeze(data_struct.data_patient.avg_vol(roi,mci_indx)));
    [f_ad,x_ad]=ksdensity(squeeze(data_struct.data_patient.avg_vol(roi,ad_indx)));
    
    figure(1), clf, hold on
    plot(x_controls, f_controls, 'g','lineWidth',1)
    plot(x_mci, f_mci, 'c','lineWidth',1)
    plot(x_ad, f_ad, 'r','lineWidth',1)
    fprintf('roi: %d %s\n', roi, data_struct.roi_name{roi})
    set(gcf,'Color',[1 1 1]);
    pause
    
end


