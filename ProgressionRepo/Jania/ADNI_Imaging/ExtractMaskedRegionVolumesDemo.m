patientDir='/Users/jaghajan/Documents/Experiments/MAPER_MASKED/Patient/';
controlDir='/Users/jaghajan/Documents/Experiments/MAPER_MASKED/Control/';
nr_roi=83; %number of regions of interest

% Extract region volumes for controls using convert3d 
[avg_vol_control vox_count_control roi_id  control_id control_id_str]=ExtractVolumesWithc3d(controlDir,nr_roi);


% Extract region volumes for patients using convert3d 
[avg_vol_patient vox_count_patient roi_id  patient_id patient_id_str]=ExtractVolumesWithc3d(patientDir,nr_roi);


%create the ipout data structure for EBDPDemonstrateADNI
data_struct.data_control.avg_vol=avg_vol_control;
data_struct.data_control.nr_vox=vox_count_control;
data_struct.data_patient.avg_vol=avg_vol_patient;
data_struct.data_patient.nr_vox=vox_count_patient;
data_struct.control_id=control_id;
data_struct.control_id_str=control_id_str;
data_struct.patient_id=patient_id;
data_struct.patient_id_str=patient_id_str;
data_struct.roi_id=roi_id;

%save the data structure
save('Mat Data/MAPER_MASKED_Data.mat','data_struct');

