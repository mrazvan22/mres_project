fileNames=dir('/cs/research/vision/camino4/camino/fonteijn/adni_freesurfer/');
count=0;
for cName=3:length(fileNames)-6
    cName
    thisName=fileNames(cName).name;
    d=importdata(['adni_freesurfer/' thisName '/stats/aseg.stats']);
    icvText=d.textdata{21};
    if isempty(strfind(icvText,'IntraCranialVol'))
        fprintf('Error no ICV found \n')
    else
        count=count+1;
        indx=strfind(icvText,',');
        icvNum=str2num(icvText(indx(3)+1:indx(4)-1));
        ICV(count).vol=icvNum;
        ICV(count).drc_study_name=thisName;
        ICV(count).adni_id=thisName(1:end-5);
        
    end
end

save('ICV.mat')