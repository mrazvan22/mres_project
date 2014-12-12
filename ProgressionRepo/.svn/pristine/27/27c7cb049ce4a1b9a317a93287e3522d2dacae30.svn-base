function SelectImagesBothMAPERandDRC

load('/Users/jaghajan/Documents/Experiments/Mat Data/MAPER_Data.mat');
load('/Users/jaghajan/Documents/Experiments/Mat Data/ICV.mat');
data_struct_old=data_struct;
clear data_struct
selected_indx_patient=[];
selected_indx_control=[];
count=0;
for i=1:length(data_struct_old.patient_id_str)
    indx=find(ismember({ICV(:).adni_id},data_struct_old.patient_id_str{i}))
    
    if ~isempty(indx)
        count=count+1;
        selected_indx_patient=[selected_indx_patient i];
        vol(count)=ICV(indx(1)).vol;
        drc_study_name{count}=ICV(indx(1)).drc_study_name;
        adni_id{count}=ICV(indx(1)).adni_id;
        
        
    end
end

data_struct.data_patient.nr_vox=data_struct_old.data_patient.nr_vox(:,selected_indx_patient);
data_struct.patient_id=data_struct_old.patient_id(selected_indx_patient);
data_struct.patient_id_str=data_struct_old.patient_id_str(selected_indx_patient);
data_struct.roi_name=data_struct_old.roi_name;
data_struct.roi_id=data_struct_old.roi_id;
data_struct.status_mat=data_struct_old.status_mat(selected_indx_patient);
data_struct.ICV_patient.vol=vol;
data_struct.ICV_patient.drc_study_name=drc_study_name;
data_struct.ICV_patient.adni_id=adni_id;

%normalize the voxel count for each region by ICV volume
aa=data_struct.data_patient.nr_vox;
bb=data_struct.ICV_patient.vol;
cc=repmat(bb,83,1);
data_struct.data_patient.avg_vol=aa./cc;





count=0;
clear vol
clear drc_study_name
clear adni_id

for i=1:length(data_struct_old.control_id_str)
    indx=find(ismember({ICV(:).adni_id},data_struct_old.control_id_str{i}))
    
    if ~isempty(indx)
        count=count+1;
        selected_indx_control=[selected_indx_control i];
        vol(count)=ICV(indx(1)).vol;
        drc_study_name{count}=ICV(indx(1)).drc_study_name;
        adni_id{count}=ICV(indx(1)).adni_id;
        
    end
end

data_struct.data_control.nr_vox=data_struct_old.data_control.nr_vox(:,selected_indx_control);
data_struct.control_id_str=data_struct_old.control_id_str(selected_indx_control);
data_struct.ICV_control.vol=vol;
data_struct.ICV_control.drc_study_name=drc_study_name;
data_struct.ICV_control.adni_id=adni_id;

aaa=data_struct.data_control.nr_vox;
bbb=data_struct.ICV_control.vol;
ccc=repmat(bbb,83,1);
data_struct.data_control.avg_vol=aaa./ccc;

save('/Users/jaghajan/Documents/Experiments/Mat Data/MAPER_ICV_NORM_DATA.mat','data_struct')


