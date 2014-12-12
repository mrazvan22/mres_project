function rt=mixTPatchVersion2Demo
 %2) function to learn a model for each patch


%SELECT DATA TYPE
%=== 1: XM2VTS ground truth data
%=== 2: Gabor filtered XM2VTS ground truth data
%=== 3: XM2VTS Light data
%=== 4: Gabor filtered XM2VTS Light data
dataType=1;

[trainingSet, imageIDTrain, name]=dataSetsTrain(dataType);

%fit gaussians to data
N_T_EST = 128;

name=[name 'Version2T' num2str(N_T_EST) '.mat'];
destinationDir='C:\PhD Research\Face Recognition\MGM\Face Mix Gauss\mixT\Patch Version\Standard Model\Mat Data\XM2VTS\';


% %libDir='C:\PhD Research\Face Recognition\MGM\Face Mix Gauss\transformations\Mat data\';
% 
% 
% 
% %============XM2VTS Data with ground truth illumination============
% %sourceDir='C:\PhD Research\Data\'
% libDir='C:\PhD Research\Code\Mosaicfaces\Standard Model\Mat data\';
% %
% % load ([libDir 'patchTrainingSets8.mat'])
% % %Load the dataIDmatrix
% % load ([sourceDir 'enhancedSets.mat']);
% 
% %============XM2VTS Data with different illumination============
% 
% sourceDir='C:\PhD Research\Code\Mosaicfaces\Appearance Model\Appearance Model on data with diff illumination\';
% destinationDir='C:\PhD Research\Code\mixT\Patch Version\Standard Model\Mat Data\';
% 
% load([sourceDir 'newLightPatchTrainingSets8.mat']);
% 
% imageIDTrain=imageID(1:size(trainPatch,2),:);
% name='Lightpatch8LibraryVersion2MixT16.mat'



%find number of total patches in each image
n_Patch=size(trainingSet,1);




%============XM2VTS Data with different illumination============

% libDir='C:\PhD Research\Code\mixT\Patch Version\Standard Model\Mat Data\Gabor\';
% sourceDir='C:\PhD Research\Code\Mosaicfaces\Appearance Model\Appearance Model on data with diff illumination\';
% load([libDir 'filteredGaborXM2VTSLightTrainingSetPatches8GrayMiddle.mat']);
% load([sourceDir 'newLightPatchTrainingSets8.mat'],'imageID');
% imageIDTrain=imageID(1:size(filteredPatchFeatures,1),:);
% 

pCount=0;
%train a model with 8 gaussians for each patch
for(cPatch=1:64)
    if(~mod(cPatch,8))
    else
        tic
        cPatch
        pCount=pCount+1;
        %data=squeeze(filteredPatchFeatures(:,cPatch,:));
        data=squeeze(trainingSet(cPatch,:,:));

        patchLibrary{pCount}=mixTDemo(data',imageIDTrain,N_T_EST );
        toc
    end

end
save([destinationDir name],'patchLibrary');