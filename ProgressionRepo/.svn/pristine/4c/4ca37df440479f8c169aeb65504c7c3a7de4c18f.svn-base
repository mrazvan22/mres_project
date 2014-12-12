function seperateControlPatientImages

%% Initialization
sourceDir='/Users/jaghajan/Documents/Experiments/MAPER_DATA/';
destinationDir='/Users/jaghajan/Documents/Experiments/MAPER_DATA_CATEGORIES/';

digitNum=4; % this is the convention used to name the IDs 4-> 0001, 0092 etc.
zeroPad{1}='';
zeroPad{2}='0';
zeroPad{3}='00';
zeroPad{4}='000';

count=0;

%% Read and process all control images
% 
% %load the controlPatientIDs
% load([sourceDir 'ControlPatientID.mat']);
% nControl=length(controlPatientID);
% 
% for cID=1:nControl
%     %extract control patient id
%     thisID=controlPatientID(cID);
%     thisIDTxt=strcat(zeroPad{digitNum-length(num2str(thisID))+1},num2str(thisID));
%     
%     %look up the relevent .nii file for this control id
%     fileNames=dir([sourceDir ['*_S_' thisIDTxt '_*']]);
%     if length(fileNames)>0
%         for cFile=1:length(fileNames)
%             copyfile([sourceDir fileNames(cFile).name],[destinationDir 'Control/' fileNames(cFile).name])
%         end
%         count=count+1;  % count how many sontrols we have the MAPER segmentation for
%         controlIDMaper(count)=thisID;
%         
%     end
%     fprintf('Procesing control %d and %d files copied \n',cID,length(fileNames))
%     
% end
% 
% %save those control patient IDs that exist in the MAPER segmentation
% save([destinationDir 'Control/ControlPatientIDMaper.mat'],'controlIDMaper');


%% Read and process all patient images

%load the controlPatientIDs
load([sourceDir 'DiagnosisInfo.mat']);
%extract patients IDs
patientIndx=find(dxcurrenBaseline~=1);
patientID=patientIDBaseline(patientIndx);
dxcurren=dxcurrenBaseline(patientIndx);
count=0;

for cID=1:length(patientID)
    %extract control patient id
    thisID=patientID(cID);
    thisDx=dxcurren(cID);
    thisIDTxt=strcat(zeroPad{digitNum-length(num2str(thisID))+1},num2str(thisID));
    
    %look up the relevent .nii file for this control id
    fileNames=dir([sourceDir ['*_S_' thisIDTxt '_*']]);
    if length(fileNames)>0
        for cFile=1:length(fileNames)
            copyfile([sourceDir fileNames(cFile).name],[destinationDir 'Patient/' fileNames(cFile).name])
        end
        count=count+1;  % count how many sontrols we have the MAPER segmentation for
        patientIDMaper(count)=thisID;
        dxcurrenMaper(count)=thisDx;
        
    end
    fprintf('Procesing control %d and %d files copied \n',cID,length(fileNames))
    
end

%save those control patient IDs that exist in the MAPER segmentation
save([destinationDir 'Patient/PatientIDMaper.mat'],'patientIDMaper','dxcurrenMaper');

