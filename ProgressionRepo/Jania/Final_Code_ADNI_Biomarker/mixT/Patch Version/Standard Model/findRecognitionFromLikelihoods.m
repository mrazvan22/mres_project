function rt=findRecognitionFromLikelihoods

%this function called faceRecognitionVersion2TDistPatchesDemo with 100 different
%gallery images and counts the number of correct recognition rate

%load appropriate patch size info
sourceDir='C:\PhD Research\Face Recognition\MGM\Face Mix Gauss\Version 3\Appearance Model on data with diff illumination\';
load ([sourceDir 'patchInfo8.mat'], 'nX','nY','nZ','patchNum');



%select data type and the library size
%SELECT DATA TYPE
%=== 1: XM2VTS ground truth data
%=== 2: Gabor filtered XM2VTS ground truth data
%=== 3: XM2VTS Light data
%=== 4: Gabor filtered XM2VTS Light data
dataType=1;
N_T_DIST=8;

%select the model type
%==1 Standard model
%== 2 Appearance Model
modelType=1; 

%load the library
[destinationDir,patchLibrary, testingSet, imageIDTest, name]=dataSetsTest(dataType,N_T_DIST,modelType);



n_Patch=length(patchLibrary);
constCov=50;
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
    [logLikelihoodModels, postModels]=faceRecognitionVersion2TDistPatchesDemoFast(k,patchLibrary,testingSet,imageIDTest,n_Patch,nX,nY,nZ,patchNum);
     
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

save ([destinationDir name]);
 
rt=n_correct;