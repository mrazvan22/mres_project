function data_struct = ...
    hd_network_Cth(dir_subj)

nr_subj = size(dir_subj, 1);
data_Cth = zeros(35, 2, size(dir_subj, 1));
data_vol = zeros(35, 2, size(dir_subj, 1));
data_subcort = zeros(49, size(dir_subj, 1));
data_ICV = zeros(size(dir_subj, 1), 1);
cnt = 1;
for pid = 1:nr_subj,
    
    file_Cth_lh = sprintf('%sstats/lh.aparc.stats', ...
        deblank(dir_subj(pid, :)));
    file_Cth_rh = sprintf('%sstats/rh.aparc.stats', ...
        deblank(dir_subj(pid, :)));
    file_subcort = sprintf('%sstats/aseg.stats', ...
        deblank(dir_subj(pid, :)));
        
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
    fid = fopen(file_subcort, 'r');
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
    C_subcort = textscan(fid, '%f%f%f%d%s%d%d%d%d%d', 49);
    fclose(fid);
    
    % Searching for ICV
    fid = fopen(file_subcort, 'r');
    flg_continue = 1;
    while flg_continue,
        
        pos_old = ftell(fid);
        sline = fgetl(fid);
        if strfind(sline, 'Intracranial Volume'),
            
            str_ICV = sline;
            flg_continue = 0;
            
        end
        
    end
    fclose(fid);
    I_comma = strfind(str_ICV, ',');
    data_ICV(pid) = str2num(str_ICV((I_comma(end-1)+1):(I_comma(end)-1)));
    
    data_vol(:, 1, cnt) = C_lh{4};
    data_vol(:, 2, cnt) = C_rh{4};
    data_Cth(:, 1, cnt) = C_lh{5};
    data_Cth(:, 2, cnt) = C_rh{5};
    data_subcort(:, cnt) = C_subcort{4};
    cnt = cnt + 1;
    
end
name_regions = C_lh{1};
name_regions_subcort = C_subcort{5};

data_struct.data_vol = data_vol;
data_struct.data_Cth = data_Cth;
data_struct.data_ICV = data_ICV;
data_struct.data_subcort = data_subcort;
data_struct.name_regions = name_regions;
data_struct.name_regions_subcort = name_regions_subcort;