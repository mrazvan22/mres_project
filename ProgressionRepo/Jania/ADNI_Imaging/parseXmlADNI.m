function parseXmlADNI

%add path to xml toolbox
addpath('/Users/jaghajan/Research Project/Matlab/xml_io_tools_2010_11_05');

%load all the .xml files for MAPER_segmentation in the main MAPER directory 
sourceDir='/Users/jaghajan/Research Project/Data/MAPER/ADNI/';
fileNames=dir([sourceDir '*MAPER_segmentation_*.xml']);
nFiles=length(fileNames);

destinationDir='/Users/jaghajan/Documents/Experiments/MAPER_DATA/';

%loop through every file
for cFile=1:nFiles
    fprintf('Processing file %d of %d\n',cFile,nFiles)
    
    %Parse the file name to extract patient ID   
    thisName=fileNames(cFile).name;
    indx=findstr(thisName,'_');
    patientID=thisName(indx(1)+1:indx(4)-1);
    
    %check if there are both Screening and Baseline images of this patient
    patientDir=dir([sourceDir patientID '/MAPER_segmentation/']);
    %if numFolders is 2 then there are both SC and BL otherwise on;y SC
    numFolders=length(find([patientDir.isdir]))-2;
    
    %if both Screening and Baseline exist
    if numFolders==2
        chosenVisit='Baseline'; % Select baseline
    else
        chosenVisit='Screening'; % Select screeninng
    end
    
    patientXmls=dir([sourceDir '*_' patientID '_MAPER_segmentation_*']);
    
    for cXml=1:length(patientXmls)
         %read the xml file into a structure array
         [thisFile thisFileDummy] = xml_read ([sourceDir patientXmls(cXml).name]);
         
         %extract the visit Identifier
         visitID=thisFile.project.subject.visit.visitIdentifier;
         
         %select if chosenVisit otherwise ignor the xml file
         if strcmp(visitID,chosenVisit)
             %extract the data acquisition date and series ID
             dataAcq=thisFile.project.subject.study.series.dateAcquired;
             seriesID=thisFile.project.subject.study.series.seriesIdentifier;
             
     
             %format filenames
             dataAcqString=strrep(dataAcq,':','_');
             dataAcqString=strrep(dataAcqString,' ','_');
             imageFile=dir([sourceDir  patientID '/MAPER_segmentation/' dataAcqString '/S' num2str(seriesID) '/*.nii']);
             imageName=imageFile.name;
             indx2=strfind(imageName,'_');
             newName=['ADNI_' patientID '_' visitID '_S' num2str(seriesID) '_' imageName(indx2(end)+1:end)];
             
             
             %store all infornation about this file
             fileInfo.patientID=patientID;
             fileInfo.visitID=visitID;
             fileInfo.dataAcq=dataAcq;
             fileInfo.seriesID=seriesID;
             fileInfo.origName=imageName;
             
             %copy the corresponding .nii image into destination folder
             copyfile([sourceDir  patientID '/MAPER_segmentation/' dataAcqString '/S' num2str(seriesID) '/' imageName],[destinationDir newName]);
             save([destinationDir strrep(newName,'.nii','_Info.mat')],'fileInfo');
         else
         end %end if chosenVisit    
         
    end % for cXml
end % for cFile