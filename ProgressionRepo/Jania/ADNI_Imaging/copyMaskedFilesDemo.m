function copyMaskedFilesDemo

%load the data structure with patient and control IDs
load('Mat Data/MAPER_DATA.mat');

%copy masked images for controls
control_id_str=data_struct.control_id_str;
copyMaskedFiles(control_id_str,'Control');

%copy masked images for patients
patient_id_str=data_struct.patient_id_str;
copyMaskedFiles(patient_id_str,'Patient');

