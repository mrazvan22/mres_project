function [data_jacc_fluid, data_jacc_ffd, data_jacc_ffd_new] = ...
    hd_network(file_data, dir_work, dir_caudate, id_roi_select)

[data_jacc_fluid, data_jacc_ffd, data_jacc_ffd_new] = ...
    extract_data(file_data, dir_work, dir_caudate, id_roi_select);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [data_jacc_fluid, data_jacc_ffd, data_jacc_ffd_new] = ...
    extract_data(file_subj, dir_work, dir_caudate, id_roi_select)

nr_roi_select = length(id_roi_select);
nr_subj = size(file_subj, 1);

data_jacc_ffd.data_mean = ...
    zeros(nr_roi_select, nr_subj, 2);
data_jacc_ffd.data_medianneg = ...
    zeros(nr_roi_select, nr_subj, 2);
data_jacc_ffd.data_medianneg_nw = ...
    zeros(nr_roi_select, nr_subj, 2);
data_jacc_ffd.data_std = ...
    zeros(nr_roi_select, nr_subj, 2);
data_jacc_ffd.data_nr_vox = ...
    zeros(nr_roi_select, nr_subj, 2);
data_jacc_fluid.data_mean = ...
    zeros(nr_roi_select, nr_subj, 2);
data_jacc_fluid.data_medianneg = ...
    zeros(nr_roi_select, nr_subj, 2);
data_jacc_fluid.data_medianneg_nw = ...
    zeros(nr_roi_select, nr_subj, 2);
data_jacc_fluid.data_std = ...
    zeros(nr_roi_select, nr_subj, 2);
data_jacc_fluid.data_nr_vox = ...
    zeros(nr_roi_select, nr_subj, 2);
data_jacc_ffd_new.data_mean = ...
    zeros(nr_roi_select, nr_subj, 2);
data_jacc_ffd_new.data_median = ...
    zeros(nr_roi_select, nr_subj, 2);

for subj = 1:nr_subj,
    
    [p, f, d1, d2] = fileparts(deblank(file_subj(subj, :)));
    I1 = strfind(f, 'patient');
    I2 = strfind(f, 'grp');
    id_patient = str2num(f((I1+7):(I2-2)));
    file_segm_nii = sprintf('%sfreesurfer/%s/mri/aparc+aseg.nii', ...
        dir_work, f);
    
    if ~exist(file_segm_nii, 'file'),

        file_segm_mgz = sprintf('%sfreesurfer/%s/mri/aparc+aseg.mgz', ...
            dir_work, f);
         unix(sprintf('mri_convert %s %s', file_segm_mgz, file_segm_nii));

    end
    if exist(file_segm_nii, 'file'),
        
        V_segm = spm_vol(file_segm_nii);
        img_segm = spm_read_vols(V_segm);
        for roi = 1:nr_roi_select

            I_roi{roi} = find(img_segm == id_roi_select(roi));

        end
        %         file_caudate = sprintf('%spt%d.nii', dir_caudate, id_patient);
        %         V_caudate = spm_vol(file_caudate);
        %         img_caudate = spm_read_vols(V_caudate);
        %         I_caudate = find(img_caudate ~= 0);

        for fu = [12 24],

            file_ffd = sprintf('%sffd_registrations/f3d_jac_pt%d_fu%d-to-baseline.nii', ...
                dir_work, id_patient, fu);
            file_fluid = sprintf('%sfluid_registrations/fluid_pt%d_fu%d-to-baseline_unsep_fixed.nii', ...
                dir_work, id_patient, fu);
            file_ffd_new = sprintf('%sF3D_BCEM_processing/jacobian/fu%d_to_bl_pt%d_nrr_jac.nii', ...
                dir_work, fu, id_patient);

            if exist(file_ffd, 'file'),
                
                V_ffd = spm_vol(file_ffd);
                img_jacc = spm_read_vols(V_ffd);
                cnt = 1;
                for roi = 1:nr_roi_select,
                    
                    I = convert_indimg(I_roi{roi}, V_segm, V_ffd);
                    if sum(img_jacc(I)),
                        
                        data_jacc_ffd.data_mean(cnt, subj, fu/12) = ...
                            mean(img_jacc(I));
                        data_jacc_ffd.data_median(cnt, subj, fu/12) = ...
                            mean(img_jacc(I));
                        data_jacc_ffd.data_std(cnt, subj, fu/12) = ...
                            std(img_jacc(I));
                        data_jacc_ffd.data_nr_vox(cnt, subj, fu/12) = ...
                            length(img_jacc(I));
                        I_neg = find(img_jacc(I) < 1);
                        data_jacc_ffd.data_medianneg(cnt, subj, fu/12) = ...
                            median(img_jacc(I(I_neg)))*(length(I_neg)/length(I));
                        data_jacc_ffd.data_medianneg_n2(cnt, subj, fu/12) = ...
                            median(img_jacc(I(I_neg)));
                        
                    end
                    cnt = cnt + 1;
                    
                end
                %                 I = convert_indimg(I_caudate, V_caudate, V_ffd);
                %                 data_jacc_ffd.data_mean_caudate(subj, fu/12) = ...
                %                     mean(img_jacc(I));
                
            end
            if exist(file_fluid, 'file'),
                
                V_fluid = spm_vol(file_fluid);
                img_jacc = exp(spm_read_vols(V_fluid));
                cnt = 1;
                for roi = 1:nr_roi_select,
                    
                    I = convert_indimg(I_roi{roi}, V_segm, V_fluid);
                    if sum(img_jacc(I)),
                        
                        data_jacc_fluid.data_mean(cnt, subj, fu/12) = ...
                            mean(img_jacc(I));
                        data_jacc_fluid.data_median(cnt, subj, fu/12) = ...
                            mean(img_jacc(I));
                        I_neg = find(img_jacc(I) < 1);
                        data_jacc_fluid.data_medianneg(cnt, subj, fu/12) = ...
                            median(img_jacc(I(I_neg)))*(length(I_neg)/length(I));
                        data_jacc_fluid.data_medianneg_nw(cnt, subj, fu/12) = ...
                            median(img_jacc(I(I_neg)));
                        data_jacc_fluid.data_std(cnt, subj, fu/12) = ...
                            std(img_jacc(I));
                        data_jacc_fluid.data_nr_vox(cnt, subj, fu/12) = ...
                            length(img_jacc(I));
                        
                    end
                    cnt = cnt + 1;
                    
                end
                %                 I = convert_indimg(I_caudate, V_caudate, V_fluid);
                %                 data_jacc_fluid.data_mean_caudate(subj, fu/12) = ...
                %                     mean(img_jacc(I));

                
            end
            if exist(file_ffd_new, 'file'),
                
                V_ffd_new = spm_vol(file_ffd_new);
                img_jacc = spm_read_vols(V_ffd_new);
                cnt = 1;
                for roi = 1:nr_roi_select,
                    
                    I = convert_indimg(I_roi{roi}, V_segm, V_ffd_new);
                    if sum(img_jacc(I)),
                        
                        data_jacc_ffd_new.data_mean(cnt, subj, fu/12) = ...
                            mean(img_jacc(I));
                        data_jacc_ffd_new.data_median(cnt, subj, fu/12) = ...
                            median(img_jacc(I));
                        data_jacc_ffd_new.data_quant(cnt, subj, fu/12) = ...
                            quantile(img_jacc(I), 0.1);
                        
                        %                         I_neg = find(img_jacc(I) < 1);
                        %                         data_jacc_fluid.data_medianneg(cnt, subj, fu/12) = ...
                        %                             median(img_jacc(I(I_neg)))*(length(I_neg)/length(I));
                        %                         data_jacc_fluid.data_medianneg_nw(cnt, subj, fu/12) = ...
                        %                             median(img_jacc(I(I_neg)));
                        %                         data_jacc_fluid.data_std(cnt, subj, fu/12) = ...
                        %                             std(img_jacc(I));
                        %                         data_jacc_fluid.data_nr_vox(cnt, subj, fu/12) = ...
                        %                             length(img_jacc(I));
                        
                    end
                    cnt = cnt + 1;
                    
                end
                %                 I = convert_indimg(I_caudate, V_caudate, V_fluid);
                %                 data_jacc_fluid.data_mean_caudate(subj, fu/12) = ...
                %                     mean(img_jacc(I));

                
            end
            
            %             V_fluid = spm_vol(file_fluid);
            %             img_jacc_fluid = spm_read_vols(V_fluid);
            %             label = roi_id(1);
            %             I = convert_indimg(I_label{label}, V_segm, V_fluid);
            %             data_caudate_fluid{subj, fu/12} = exp(img_jacc_fluid(I));
            %             V_ffd = spm_vol(file_ffd);
            %             img_jacc_ffd = spm_read_vols(V_ffd);
            %             label = roi_id(1);
            %             I = convert_indimg(I_label{label}, V_segm, V_ffd);
            %             data_caudate_ffd{subj, fu/12} = img_jacc_ffd(I);

            fprintf('fu: %d\tsubject: %d\n', fu, subj)

        end
        
    end
    
end
    