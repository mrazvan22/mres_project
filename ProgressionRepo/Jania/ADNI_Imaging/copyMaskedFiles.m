function copyMaskedFiles(patientID,type)

sourceDir='/Users/jaghajan/Research Project/Data/MAPER/ADNI/';
infoDir=['/Users/jaghajan/Documents/Experiments/MAPER_DATA_CATEGORIES/' type '/'];
destinationDir=['/Users/jaghajan/Documents/Experiments/MAPER_MASKED/' type '/'];

% for each patient/control copy the corresponding mask files
for cID=1:length(patientID)
    cID
    %find and load the corresonding control info file
    fileName=dir([infoDir '*' patientID{cID} '*.mat']);
    load([infoDir fileName.name]);
    %extract data aquisiton time and series ID to select the right folder
    dataAcq=fileInfo.dataAcq;
    dataAcqString=strrep(dataAcq,':','_');
    dataAcqString=strrep(dataAcqString,' ','_');
    seriesIDString=['S' num2str(fileInfo.seriesID)];
             
    copyfile([sourceDir patientID{cID} '/MAPER_segmentation,_masked/' dataAcqString '/' seriesIDString '/*.nii'],  destinationDir);
    
end

