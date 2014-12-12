function rt=findRecognitionFromLikelihoodsTest

%this function called faceRecognitionVersion2TDistPatchesDemo with 100 different
%gallery images and counts the number of correct recognition rate

destinationDir='C:\PhD Research\Code\mixT\Patch Version\Standard Model\Mat Data\';

%name='LightPatch8Version2MixT32';
 
%load the library and testing set 
load([destinationDir 'Lightpatch8LibraryVersion2MixT8.mat'],'patchLibrary');
%load([destinationDir 'filteredGaborXM2VTSLightTestingSetPatches8GrayMiddleIllumCorrect.mat']);
 

% %=========XM2VTS data with ground truth illumination ==================
% libDir='C:\PhD Research\Code\mixT\Patch Version\Standard Model\Mat Data\'; 
 sourceDir='C:\PhD Research\Code\Mosaicfaces\Standard Model\Mat data\';
 sourceDir2='C:\PhD Research\Data\';
% % %load the library
%  load([libDir 'patch8LibraryVersion2MixT8.mat'],'patchLibrary');
% % 
% %load the test set and the probes
 load ([sourceDir 'patchTestingSets8.mat']);
% 
 %Load the dataIDmatrix
 load ([sourceDir2 'enhancedSets.mat']);

% % % % %=========XM2VTS data with different illumination ==================
%     sourceDir='C:\PhD Research\Code\Mosaicfaces\Appearance Model\Appearance Model on data with diff illumination\';
% % % % %  
%     load([sourceDir 'newLightPatchTestingSets8.mat']);
% % 
%     name='Lightpatch8Version2MixT8.mat';
% %   

 
 
 
%find number of total patches in each image
%n_Patch=size(testPatch,1);
%n_Patch=size(filteredPatchFeatures,2);
%load the size of the image and number of the patches 
 load ([sourceDir 'patchInfo8.mat'], 'nX','nY','nZ','patchNum');
%imageIDTest=imageID(781:end,:);
%imageIDTest=dataIDMatrix;

n_Patch=length(patchLibrary);
constCov=100;
constCov_flag=0;
 
if(constCov_flag)
    for(cPatch=1:length(patchLibrary))
        l=patchLibrary{cPatch};
        [N_IND N_APP]=size(l);
            for(cInd=1:N_IND)
                for(cApp=1:N_APP)
                    l(cInd,cApp).cov=l(cInd,cApp).cov*constCov;
                end
            end
        patchLibrary{cPatch}=l;    
    end
end

n_correct=0;
incorrect=[];



for(k=1:100)
    tic
    k
    [logLikelihoodModels, postModels]=faceRecognitionVersion2TDistPatchesDemoFastTest(k,patchLibrary,testPatch,imageIDTest,n_Patch,nX,nY,nZ,patchNum);
     
    %[gallery, probe]=showData(k,patchLibrary,testPatch,imageIDTest,n_Patch,nX,nY,nZ,patchNum);
    
    likelihoods(k,:,:)=logLikelihoodModels;
    
    
    s=sum(logLikelihoodModels);
     s2=find(s==max(s));
     est_indx(k)=s2(1);
     
     posterior=s-(max(s));
     posterior=exp(posterior);
     posterior=posterior/sum(posterior);
     
     if(est_indx(k)==k)
         n_correct=n_correct+1;
     else
         incorrect=[incorrect k];
     end
     toc
end

%save ([destinationDir name]);
 
rt=n_correct;