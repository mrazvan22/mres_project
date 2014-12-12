%% Step by step Image and data extraction from MAPER ADNI


%% Process the MAPER Segmented images

%1) Parse the Xml files in the MAPER ADNI and copy all the baseline images in
%a seperate directory > /Users/jaghajan/Documents/Experiments/MAPER_DATA
parseXmlADNI;

%2) Parse the CSV turned Xls file to store the baseline and control patient ID
%and Diagnosis info. 
parseXlsADNI;

%3) Run the shell script (from Terminal) to run convert3d command line tool to estimate
%voxel count for each region
run3d.sh

%4) Separate the control and patient Images and store them in directory > 
% /Users/jaghajan/Documents/Experiments/MAPER_DATA_CATEGORIES
seperateControlPatientImages;

%5)Extract Average Volume measures for each region for control and patient
ExtractRegionVolumesDemo;


%% Process the MAPER Masked images
%6) Copy the masked images of patient/control to a separate directory >
%/Users/jaghajan/Documents/Experiments/MAPER_MASKED/
copyMaskedFilesDemo;
