function data_struct = ...
    ad_network_Cth(dir_control, dir_patient)

work_dir = '/cs/research/vision/camino1/camino/fonteijn/forDannyAndHubert';

file_lut = [work_dir '/FreeSurferColorLUT.txt'];
[label_id, label_name, d1,d2, d3, d4] = ...
    textread(file_lut, '%d%s%d%d%d%d');
nr_labels = length(label_id);

[data_Cth_control, data_vol_control, name_regions] = ...
    extract_data(dir_control, 'c', label_id, work_dir);
[data_Cth_patient, data_vol_patient] = ...
    extract_data(dir_patient, 'p', label_id, work_dir);

data_struct.data_Cth_control = data_Cth_control;
data_struct.data_Cth_patient = data_Cth_patient;
data_struct.data_vol_control = data_vol_control;
data_struct.data_vol_patient = data_vol_patient;
data_struct.name_regions = name_regions;
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [data_Cth, data_vol, name_regions] = ...
    extract_data(dir_subj, str_corp, label_id, work_dir)

for tp = 1:size(dir_subj, 1),
    
    subj_id_str(tp, :) = dir_subj(tp, (end-6):(end-1));
    subj_id(tp) = str2num(subj_id_str(tp, 2:3));
    subj_tp(tp) = str2num(subj_id_str(tp, 5:6));
    
end

subj_id_unique = unique(subj_id);
data_Cth = zeros(35, 2, size(dir_subj, 1) - length(subj_id_unique));
data_vol = zeros(35, 2, size(dir_subj, 1) - length(subj_id_unique));
cnt = 1;
for pid = subj_id_unique,
    
    I_subj = find(subj_id == pid);
    subj_tp_local = subj_tp(I_subj);
    for t = 1:max(subj_tp_local),
                
        file_Cth_lh = sprintf('%s/freesurfer/%s%.2dt%.2d/stats/lh.aparc.stats', ...
            work_dir, str_corp, pid, t);
        file_Cth_rh = sprintf('%s/freesurfer/%s%.2dt%.2d/stats/rh.aparc.stats', ...
            work_dir, str_corp, pid, t);
        
        fid = fopen(file_Cth_lh, 'r');
        flg_continue = 1;
        while flg_continue,
            
            pos_old = ftell(fid);
            sline = fgetl(fid);
            if sline(1) ~= '#',
                
                pos = pos_old;
                flg_continue = 0;
                
            end
            
        end
        fseek(fid, pos, 'bof');
        C_lh = textscan(fid, '%s%d%d%d%f%f%f%f%d%f', 35);
        fclose(fid);
        fid = fopen(file_Cth_rh, 'r');
        flg_continue = 1;
        while flg_continue,
            
            pos_old = ftell(fid);
            sline = fgetl(fid);
            if sline(1) ~= '#',
                
                pos = pos_old;
                flg_continue = 0;
                
            end
            
        end
        fseek(fid, pos, 'bof');
        C_rh = textscan(fid, '%s%d%d%d%f%f%f%f%d%f', 35);
        fclose(fid);
       
        data_vol(:, 1, cnt) = C_lh{4};
        data_vol(:, 2, cnt) = C_rh{4};
        data_Cth(:, 1, cnt) = C_lh{5};
        data_Cth(:, 2, cnt) = C_rh{5};
        cnt = cnt + 1;
        
    end
    
end
name_regions = C_lh{1};