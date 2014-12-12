function rt=mixGaussPatchVersion3Demo

%2) function to learn a model for each patch 
sourceDir='C:\PhD Research\Code\Mosaicfaces\Appearance Model\Appearance Model on data with diff illumination\';

destinationDir='C:\PhD Research\Code\Mosaicfaces\Appearance Model\Appearance Model on data with diff illumination\';


N_IND=25;
N_APP=4;
%sourceDir='C:\PhD Research\Face Recognition\MGM\Face Mix Gauss\Version 3\Mat Data\FERRET\';
%load the imageIdtrain
% load([sourceDir 'Pix_FERET_TrainSet.mat'])
% imageID=dataIDMatrix;         
%Load patch training set dataIDMatrix
load([sourceDir 'newLightPatchTrainingSets8.mat']);
% load([sourceDir 'filteredGaborFERETTrainingSetPatches16GrayMiddle.mat']);
%find number of total patches in each image
%n_Patch=size(filteredPatchFeatures,2);
%find number of total patches in each image
n_Patch=size(trainPatch,1);

imageID=imageID(1:size(trainPatch,2),:);
%train a model with 8 gaussians for each patch
for(cPatch=1:40)
    tic
    if (~mod(cPatch,8))
    else
        cPatch
        data=squeeze(trainPatch(cPatch,:,:));
        %data=squeeze(filteredPatchFeatures(:,cPatch,:));
        patchLibrary{cPatch}=mixTPatchBasedVersion3(data,imageID,N_IND,N_APP);
    end
    toc
end
save([destinationDir 'XM2VTSDiffIllumPatch8LibraryVersion3T254Fast.mat'],'patchLibrary');