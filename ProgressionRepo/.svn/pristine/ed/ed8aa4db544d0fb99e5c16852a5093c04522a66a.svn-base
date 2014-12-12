%step by step FreeSurfer atrophy extraction

%copy the relevant freesurfer ,mgz files to comic100
$ scp -r /Users/jaghajan/Research_Project/Data/FreeSurfer/Longitudinal/ jaghajan@comic100:/home/jaghajan/FreeSurferImages/Longitudinal

%write a scipt to convert all mgz files in convertMGZ.sh Do for both
%patient and control
$ chmod =x convertMGZ.sh
$ sh convertMGZ.sh

%copy the .nii.gz files back on local machine run from destination dir
$ scp jaghajan@comic100:/home/jaghajan/FreeSurferImages/Longitudinal/Patient/*.nii.gz .
$ scp jaghajan@comic100:/home/jaghajan/FreeSurferImages/Longitudinal/Contol/*.nii.gz .

%unzip all the .nii.gz files to .nii from both control and patient
%directories
$ gunzip *.gz

%run convert3d to get the resgion stats from the scripts run3dFreeSurfer.sh
sh run3dFreeSurfer.sh  %for both control and patient