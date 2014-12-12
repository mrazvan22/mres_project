%% ===================== load the data ==============================
close all

file_data = '/Users/jaghajan/Documents/Experiments/Mat Data/MAPER_ICV_NORM_DATA.mat'
load(file_data);

% divide the data into 4 subgroups
subgroup_indx_patinet{1}=[1:119];
subgroup_indx_patinet{2}=[120:238];
subgroup_indx_patinet{3}=[239:357];
subgroup_indx_patinet{4}=[358:478];

subgroup_indx_control{1}=[1:49];
subgroup_indx_control{2}=[50:98];
subgroup_indx_control{3}=[99:147];
subgroup_indx_control{4}=[148:197];

version_likelihood=6;
n_repeat=10;
div_size=[0.1,0.3,0.5]

% figure
% count=0;
% for i_var=1:2:83
%     count=count+1;
%     [hist_c_all, x_c_all] = ksdensity(data_struct.data_control.avg_vol(i_var,:));
%     [hist_p_all, x_p_all] = ksdensity(data_struct.data_patient.avg_vol(i_var,:));
%     
%     subplot(6,7,count), hold on
%     plot(x_c_all, hist_c_all, 'g');
%     plot(x_p_all, hist_p_all, 'r'); 
%     title(data_struct.roi_name{i_var})
%     set(gcf,'Color',[1 1 1]);
% end

for cGroup=1:1
    cGroup
    this_group_data_struct.data_patient.nr_vox=data_struct.data_patient.nr_vox(:,subgroup_indx_patinet{cGroup});
    this_group_data_struct.data_patient.avg_vol=data_struct.data_patient.avg_vol(:,subgroup_indx_patinet{cGroup});
    this_group_data_struct.patient_id=data_struct.patient_id(subgroup_indx_patinet{cGroup});
    this_group_data_struct.patient_id_str={data_struct.patient_id_str{subgroup_indx_patinet{cGroup}}};
    this_group_data_struct.roi_name=data_struct.roi_name;
    this_group_data_struct.roi_id=data_struct.roi_id;
    this_group_data_struct.status_mat=data_struct.status_mat(subgroup_indx_patinet{cGroup});
    this_group_data_struct.ICV_patient.vol=data_struct.ICV_patient.vol(subgroup_indx_patinet{cGroup});
    this_group_data_struct.ICV_patient.drc_study_name={data_struct.ICV_patient.drc_study_name{subgroup_indx_patinet{cGroup}}};
    this_group_data_struct.ICV_patient.adni_id={data_struct.ICV_patient.adni_id{subgroup_indx_patinet{cGroup}}};
    
    this_group_data_struct.data_control.nr_vox=data_struct.data_control.nr_vox(:,subgroup_indx_control{cGroup});
    this_group_data_struct.data_control.avg_vol=data_struct.data_control.avg_vol(:,subgroup_indx_control{cGroup});
    this_group_data_struct.control_id_str={data_struct.control_id_str{subgroup_indx_control{cGroup}}};
    this_group_data_struct.ICV_control.vol=data_struct.ICV_control.vol(subgroup_indx_control{cGroup});
    this_group_data_struct.ICV_control.drc_study_name={data_struct.ICV_control.drc_study_name{subgroup_indx_control{cGroup}}};
    this_group_data_struct.ICV_control.adni_id={data_struct.ICV_control.adni_id{subgroup_indx_control{cGroup}}};

    for cSize=1:length(div_size)
        cSize
    robustnessTest2ImagingSubGroups(div_size(cSize),version_likelihood,n_repeat,this_group_data_struct,cGroup)
    end
end