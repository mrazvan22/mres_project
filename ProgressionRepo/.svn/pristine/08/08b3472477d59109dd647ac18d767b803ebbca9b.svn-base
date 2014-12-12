% Bit of code that converts atrophy orderings into images
% For this I'm using the segmentation of c01t01
% I need three images here: the anot+aseg image from freesurfer and the
% manually segmented left and right hippocampus. The last two you can find
% in
% /cs/research/vision/camino1/camino/fonteijn/forDannyAndHubert/data
% The ordering is: 1 = aseg+annot, 2 = leftHippocampus, 3 =
% rightHippocampus

roi_select = [17 1001:1034 53 2001:2034]; % Freesurfer indicators that are hippocampus (17 = left, 53 = right) and the cortical areas
nr_roi_select = length(roi_select);

file_segm = spm_select(3, 'image');
V_fs_parc = spm_vol(file_segm(1, :));
img_fs_parc = spm_read_vols(V_fs_parc);
V_hipp = spm_vol(file_segm(2:3, :));
img_hipp = spm_read_vols(V_hipp);
% Deleting the hippocampus from the original freesurfer segmentation
img_fs_parc(img_fs_parc == 17) = 0;
img_fs_parc(img_fs_parc == 53) = 0;
I_hipp_lh = find(img_hipp(:, :, :, 1) > 0);
% convert_indimg is a short wrapper I road to convert coordinates from
% images in different coordinate frames to each other.
I_hipp_lh_parc = convert_indimg(I_hipp_lh, V_hipp(1), V_fs_parc);
I_hipp_rh = find(img_hipp(:, :, :, 2) > 0);
I_hipp_rh_parc = convert_indimg(I_hipp_rh, V_hipp(2), V_fs_parc);
img_fs_parc(I_hipp_lh_parc) = 17;
img_fs_parc(I_hipp_rh_parc) = 53;

avPos_control = load('ControlsBothHemiSeqRegionsOutput');
avPos_patient = load('PatientsBothHemiSeqRegionsOutput');

img_order_control = zeros(size(img_fs_parc));
img_order_patient = zeros(size(img_fs_parc));
img_segm = zeros(size(img_fs_parc));
% Now going through all the ROIs and filling in the temporal value
for roi = 1:nr_roi_select
    
    I_roi = find(img_fs_parc == roi_select(roi));
    img_segm(I_roi) = roi;
    img_order_control(I_roi) = avPos_control.avPos(roi);
    img_order_patient(I_roi) = avPos_patient.avPos(roi);
    
end

V_order_control = V_fs_parc;
V_order_patient = V_fs_parc;
V_order_control.fname = 'img_ordering_control.nii';
V_order_patient.fname = 'img_ordering_patient.nii';
spm_create_vol(V_order_control);
spm_create_vol(V_order_patient);
spm_write_vol(V_order_control, img_order_control);
spm_write_vol(V_order_patient, img_order_patient);