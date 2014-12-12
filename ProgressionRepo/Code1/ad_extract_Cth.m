clear all
close all

work_dir = '/cs/research/vision/camino1/camino/fonteijn/forDannyAndHubert';
dir_control = spm_select([1 100], 'dir', 'Give controls');
dir_patient = spm_select([1 100], 'dir', 'Give patients');

for tp = 1:size(dir_patient, 1),
    
    patient_id_str(tp, :) = dir_patient(tp, (end-6):(end-1));
    patient_id(tp) = str2num(patient_id_str(tp, 2:3));
    patient_tp(tp) = str2num(patient_id_str(tp, 5:6));
    
end
for tp = 1:size(dir_control, 1),
    
    control_id_str(tp, :) = dir_control(tp, (end-6):(end-1));
    control_id(tp) = str2num(control_id_str(tp, 2:3));
    control_tp(tp) = str2num(control_id_str(tp, 5:6));
    
end

% Removing all first time points from the list.

[data_struct] = ...
    ad_network_Cth(dir_control, dir_patient);
data_struct.patient_id = patient_id;
data_struct.patient_tp = patient_tp;
data_struct.control_id = control_id;
data_struct.control_tp = control_tp;
save data_AD_CthVol data_struct

% data_Cth_control = squeeze(mean(data_struct.data_Cth_control, 2));
% data_Cth_patient = squeeze(mean(data_struct.data_Cth_patient, 2));
% data_vol_control = squeeze(mean(data_struct.data_vol_control, 2));
% data_vol_patient = squeeze(mean(data_struct.data_vol_patient, 2));
% 
% cnt1 = 1;
% cnt2 = 1;
% for pid = unique(patient_id),
%     
%     I_subj = find(patient_id == pid);
%     Cth_baseline = data_Cth_patient(:, I_subj(patient_tp(I_subj) == 1));
%     for tp = 2:max(patient_tp(I_subj)),
%         
%         Cth_t = data_Cth_patient(:, I_subj(patient_tp(I_subj) == tp));
%         data_Cth_patient_t2baseline(:, cnt1) = ...
%             (Cth_t)./Cth_baseline;
%         cnt1 = cnt1 + 1;
%         
%     end
%     vol_baseline = data_vol_patient(:, I_subj(patient_tp(I_subj) == 1));
%     for tp = 2:max(patient_tp(I_subj)),
% 
%         vol_t = data_vol_patient(:, I_subj(patient_tp(I_subj) == tp));
%         data_vol_patient_t2baseline(:, cnt2) = ...
%             (vol_t)./vol_baseline;
%         cnt2 = cnt2 + 1;
% 
%     end
%     
% end
% 
% cnt1 = 1;
% cnt2 = 1;
% for pid = unique(control_id),
%     
%     I_subj = find(control_id == pid);
%     Cth_baseline = data_Cth_control(:, I_subj(control_tp(I_subj) == 1));
%     for tp = 2:max(control_tp(I_subj)),
%         
%         Cth_t = data_Cth_control(:, I_subj(control_tp(I_subj) == tp));
%         data_Cth_control_t2baseline(:, cnt1) = ...
%             (Cth_t)./Cth_baseline;
%         cnt1 = cnt1 + 1;
%         
%     end
%     vol_baseline = data_vol_control(:, I_subj(control_tp(I_subj) == 1));
%     for tp = 2:max(control_tp(I_subj)),
% 
%         vol_t = data_vol_control(:, I_subj(control_tp(I_subj) == tp));
%         data_vol_control_t2baseline(:, cnt2) = ...
%             (vol_t)./vol_baseline;
%         cnt2 = cnt2 + 1;
% 
%     end
% 
% end
% 
% atrophy_patient = data_vol_patient_t2baseline';
% atrophy_control = data_vol_control_t2baseline';
% [numMeas numEvents] = size(atrophy_patient);
% p1 = PrecedenceMatrix1(atrophy_patient);
% p2 = PrecedenceMatrix2(atrophy_patient, atrophy_control);
% 
% p = p2;
% %% Initialize the model
% 
% % Construct starting path.
% [junk samplePath] = sort(rand(1,numEvents));
% 
% % Evaluate likelihood
% pathlik = 0;
% for i=1:(numEvents-1)
%     for j=(i+1):numEvents
%         pathlik = pathlik + log(p(samplePath(i), samplePath(j)));
%     end
% end
% 
% %% Run the MCMC
% curPath = samplePath;
% curlik = pathlik;
% 
% burnin = 5000;
% interval = 50;
% numSamples = 2000;
% iterations = interval*numSamples+burnin;
% samples = zeros(numSamples, numEvents);
% liks = zeros(numSamples, 1);
% 
% % Run MCMC
% sno = 1;
% for it=1:iterations
%     
%     % Perturb current path
%     ind1 = fix(rand(1,1)*numEvents)+1;
%     ind2 = ind1;
%     while(ind2==ind1)
%         ind2 = fix(rand(1,1)*numEvents)+1;
%     end
%     
%     % Switch entries
%     newPath = curPath;
%     newPath(ind1) = curPath(ind2);
%     newPath(ind2) = curPath(ind1);
%     
%     % Evaluate likelihood
%     newlik = 0;
%     for i=1:(numEvents-1)
%         for j=(i+1):numEvents
%             newlik = newlik + log(p(newPath(i), newPath(j)));
%         end
%     end
% 
%     if(newlik>curlik)
%         acc = 1;
%     else
%         a = exp(newlik - curlik);
%         acc = rand<a;
%     end
%     
%     if(acc)
%         curlik = newlik;
%         curPath = newPath;
%     end
% 
%     if(it>burnin && mod(it-burnin, interval)==0)
%         samples(sno,:) = curPath;
%         liks(sno) = curlik;
%         sno = sno+1;
%         display(sprintf('It: %i. Lik: %f', it, curlik));
%     end
% end
% 
% [uliks ulikinds n] = unique(liks);
% samples(ulikinds,:)
% sumord=zeros(1,numEvents);
% sumprobs = 0;
% for i=1:length(ulikinds)
%     [a b] = sort(samples(ulikinds(i),:));
%     prob = exp(uliks(i));
%     sumord = sumord+b*prob;
%     sumprobs = sumprobs + prob;
% end
% avPos = sumord/sumprobs
% 
% % Making an image of the results...
% file_segm = sprintf('%smri/aparc+aseg.nii', ...
%     deblank(dir_patient(1, :)));
% V_segm = spm_vol(file_segm);
% img_segm = spm_read_vols(V_segm);
% 
% % FreesurferColorLUT contains the Freesurfer segment names and their values 
% file_lut = [work_dir '/FreeSurferColorLUT.txt'];
% [label_id, label_name, d1,d2, d3, d4] = ...
%     textread(file_lut, '%d%s%d%d%d%d');
% 
% img_order = zeros(V_segm.dim);
% for lhrh = 0,
%     
%     for region1 = 1:length(name_regions)
%         
%         if lhrh == 0,
%             
%             name_local = ['lh-' name_regions{region1}];
%             
%         elseif lhrh == 1,
%             
%             name_local = ['rh-' name_regions{region1}];
%             
%         end
%         flg_continue = 1;
%         region2 = 1;
%         while (flg_continue && (region2 <= length(label_name))),
%             
%             if strfind(label_name{region2}, name_local),
%                 
%                 label_id_local = label_id(region2);
%                 I_region = find(img_segm == label_id_local);
%                 img_order(I_region) = ...
%                     avPos(region1 + lhrh*length(name_regions));
%                 flg_continue = 0;
%                 
%             end
%             region2 = region2 + 1;
%             
%         end
%         
%     end
%     
% end
% V_order = V_segm;
% V_order.fname = 'img_order_patient_vol_precmat2.nii';
% spm_create_vol(V_order);
% spm_write_vol(V_order, img_order);
