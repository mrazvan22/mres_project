%a function to perform face recognition on 100 testing data on version 2
%model i.e. only identity hidden variables.
function [pr_models, postModel]=faceRecognitionVersion2TDistPatchesDemo(k,patchLibrary,testPatch,imageIDTest,n_Patch,nX,nY,nZ,patchNum)




%get the dimensions of each patch
patch_sizeX=nX/patchNum;
patch_sizeY=nY/patchNum;

draw_flag1=0;
draw_flag2=0;
pr_models=zeros(n_Patch,100);

pCount=0;
%select patch region and do model selection for each patch region
for(cPatch=1:40)
     
    if (~ mod(cPatch,8)) 
    else
    
         pCount=pCount+1;
        %collect the testing set for
        testingSet=squeeze(testPatch(:,cPatch,:));
        
        gaussEst=patchLibrary{cPatch};
        
        %gaussEst=patchLibrary{pCount};
        
        
         
        %select testing set
        for(cPerson=1:size(imageIDTest,2))
            indices=find(imageIDTest(:,cPerson));
            gallery(:,cPerson)=testingSet(indices(1),:);
            probe(:,cPerson)=testingSet(indices(2),:);
        end

        %initialization of variables
        n_model=size(gallery,2);
        n_data=size(gallery,2);
        postModel=ones(1,n_model);

        %save ([libDir 'galleryProbe2Patch' num2str(cPatch) '.mat'],'gallery','probe');
        %load ([testDir 'galleryProbe.mat'],'gallery','probe');

        n_correct=0;
        count=1;
         
%         tic
%         %find all single terms
%         for (m=1:n_model)
%             S(m)=getSingleTermTDist(gallery(:,m),gaussEst);
%         end
%         toc
%         iii=3
%         tic


        %=============== Find single term=======================
        n_hidden=size(gaussEst,2);
        for(j=1:n_hidden)
            term(j,:)=getLogTProb(gallery,gaussEst(j).mean,gaussEst(j).cov,gaussEst(j).dof)+log(gaussEst(j).prior);
        end

        %rescale log probabilities by their maximim
        maxTerm=max(term(:));
        term=term-maxTerm;

        %do summation
        summedTerm=sum(exp(term));
        summedTerm(summedTerm==0)=realmin;
        %convert back to log probabilities
        term=log(summedTerm)+maxTerm;
        S=term;  
        %======================================================
        
        
        %select the probe image
        XP=probe(:,k);

        if(draw_flag1)
            %find the most likeliy patch in the library that explain
            %the probe and the gallery image
            [mostLikelyLibPatch dataLogLike]=findMostLikelyLibPatch(XP,cPatch);
            reconsPatch{cPatch}=gaussEst(mostLikelyLibPatch).mean;
            %plot the gallery and probe and the most likely library patch
            subplot(2,2,1);imshow(reshape(uint8(gallery(:,k)),patch_sizeX,patch_sizeY,nZ));
            subplot(2,2,2); imshow(reshape(uint8(XP),patch_sizeX,patch_sizeY,nZ));
            subplot(2,2,3:4); imshow(reshape(uint8(gaussEst(mostLikelyLibPatch).mean),patch_sizeX,patch_sizeY,nZ));
            drawnow;
        end


%         tic
%         %calculate the likelihodd of the models given the probe
%         for(i=1:n_model)
%             pr_models(pCount,i)=getProbModelTDist(gallery, XP,gaussEst,S,i);
%         end
%         toc
        
        %=======================Find double term============================

        
        for(j=1:n_hidden)
            term(j,:)=getLogTProb(gallery,gaussEst(j).mean,gaussEst(j).cov,gaussEst(j).dof)+ ...
                getLogTProb(XP,gaussEst(j).mean,gaussEst(j).cov,gaussEst(j).dof)+log(gaussEst(j).prior);
        end

        
         %rescale log probabilities by their maximim
        maxTerm=max(term(:)); 
        term=term-maxTerm;

        %do summation
        summedTerm=sum(exp(term));
        summedTerm(summedTerm==0)=realmin;
        
        %convert back to log probabilities
        term=log(summedTerm)+maxTerm;
        
        for(i=1:n_model)
            pr_models(pCount,i)=(sum(S)+term(i))-S(i);
        end
        
         
        
        %==================================================================
        
        %find the posterior probability of the model as the product of
        %thepoeteriors of all patches
        postModel=postModel+pr_models(cPatch,:);
        if(draw_flag2)
            %plot the likelihood of the models given this pathces
            subplot(3,2,cPatch);plot(pr_models(cPatch,:),'r-');
            hold on;
            drawnow;
        end
        %end
    end
     
end
% [X,Y]=meshgrid([1:n_Patch],[1:n_model]);
% mesh(X,Y,pr_models);
% plot(pr_models)

%Rescaling the probabilities
maxPrm=max(postModel(:));
postModel=postModel-maxPrm;

postModel=exp(postModel);

posterior_models=postModel/sum(postModel);
%Finding the estimated identity/class/index


%save([libDir 'recognitionRatePatch2.mat'],'n_correct','postModel');

if(draw_flag2)
    %plot the posterior probabilities of all models
    hold on;
    subplot(3,2,5:6);plot(posterior_models,'r-');
end
%draw the reconstructed image
if(draw_flag1)
    sourceDir='C:\PhD Research\Face Recognition\Probabilistic FR\Enhanced Images\';
    load ([sourceDir 'enhancedSets.mat']);

    %plot reconstructed face image
    figure;
    reconX=buildImage(patchNum,reconsPatch,patch_sizeX,patch_sizeY,nZ);
    indices=find(imageIDTest(:,k));
    faceXp=testingSet(indices(2),:);
    faceG=testingSet(indices(1),:);
    subplot(2,2,1); imshow(reshape(uint8(faceG),70,70,3));
    title('Gallery Image');
    hold on;
    subplot(2,2,2); imshow(reshape(uint8(faceXp),70,70,3));
    title('Probe Image');
    subplot(2,2,3:4); imshow(uint8(reconX));
    title('Reconstrcuted Image');
    drawnow;
end


        
;