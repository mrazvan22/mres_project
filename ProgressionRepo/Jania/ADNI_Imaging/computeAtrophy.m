function computeAtrophy(sourceDir,nr_roi)

%select all screening images
fileNames=dir([sourceDir, '*Scr*.txt']);
scrVoxCount=zeros(nr_roi,length(fileNames));
atrophy=zeros(nr_roi,length(fileNames));

for cName=1:length(fileNames)
    %load the screening file
    thisName=fileNames(cName).name;
    thisImStats=importdata([sourceDir thisName]);
    
    %extract the patient/control ID
    %patientID=
    indx=strfind(thisName,'_');
    patientID(cName)=str2num(thisName(indx(2)+1:indx(3)-1));
    patientIDFull{cName}=thisName(1:indx(3)-1);
    
    %extract screening image voxel count
    scrVoxCount=thisImStats.data(:,6);
    
    %extract m12 region voxel count and compute absolute atrophy
    m12Imstats=importdata([sourceDir strrep(thisName,'Scr','m12')]);
    atropy(:,cName)=abs(m12Imstats.data(:,6)-scrVoxCount);
end
    

