function [regionVolume regVoxCount regionID patientID patientIDFull]=ExtractVolumesWithc3d(sourceDir,nr_roi)

%% Initialization
fileNames=dir([sourceDir '*.txt']);
regionVolume=zeros(nr_roi,length(fileNames));
regVoxCount=zeros(nr_roi,length(fileNames));

%% Extract and compute region volumes
%load each txt file and extract the regionID i.e. and the pixel count in the first and last columns

for cName=1:length(fileNames)
    cName
    thisName=fileNames(cName).name;
    thisImgStats=importdata([sourceDir thisName]);
    
    %extract the patient/control ID
    %patientID=
    indx=strfind(thisName,'_');
    patientID(cName)=str2num(thisName(indx(3)+1:indx(4)-1));
    patientIDFull{cName}=thisName(indx(1)+1:indx(4)-1);
    
    
    %extract the regionID Note: exclude the first entry i.e. regionID=0 is the
    %background region
    regionID=thisImgStats.data(2:end,1);
    
    %extract the pixel count for each region Note exclude the first entry is the
    %pixel count for the background region
    regionPixCount=thisImgStats.data(2:end,6);
    
    % find out the average volume for each region averaged by brain size
    brainSize=sum(regionPixCount);
    
    regionVolume(:,cName)=regionPixCount/brainSize;
    regVoxCount(:,cName)=regionPixCount;
end