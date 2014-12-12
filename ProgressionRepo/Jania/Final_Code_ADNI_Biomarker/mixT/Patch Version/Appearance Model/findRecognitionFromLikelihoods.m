function rt=findRecognitionFromLikelihoods

%this function called faceRecognitionVersion2PatchesDemo with 100 different
%gallery images and counts the number of correct recognition rate

destinationDir='C:\PhD Research\Code\Mosaicfaces\Appearance Model\Appearance Model on data with diff illumination\'; 
name='XM2VTSDiffIllumPatch8Version3GaussT254.mat';

libDir='C:\PhD Research\Code\Mosaicfaces\Appearance Model\Appearance Model on data with diff illumination\'; 
%load([libDir 'newLightPatch8LibraryVersion3G82.mat'],'patchLibrary');
load([libDir 'XM2VTSDiffIllumPatch8LibraryVersion3T254Fast.mat'],'patchLibrary');

constCov=30;
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

N_IND=8;
N_APP=2;

n_correct=0;
incorrect=[];
for(k=1:100)
    tic
    k
    logLikelihoodModels=faceRecognitionVersion3PatchesTDemo(k,patchLibrary);
    
    logLikeAllProbe(k,:,:)=logLikelihoodModels;
    s=sum(logLikelihoodModels);
     est_indx(k)=find(s==max(s));
     
%      Visualization(k,est_indx);
%      drawnow
      if(est_indx(k)==k)
         n_correct=n_correct+1;
     else
         incorrect=[incorrect k];
      end
     toc
end

save ([destinationDir name]);
 
rt=n_correct;