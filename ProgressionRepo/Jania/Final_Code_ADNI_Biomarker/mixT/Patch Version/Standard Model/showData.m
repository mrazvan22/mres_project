%a function to perform face recognition on 100 testing data on version 2
%model i.e. only identity hidden variables.
function [gallery, probe]=showData(k,patchLibrary,testPatch,imageIDTest,n_Patch,nX,nY,nZ,patchNum)




%get the dimensions of each patch
patch_sizeX=nX/patchNum;
patch_sizeY=nY/patchNum;

draw_flag1=0;
draw_flag2=0;
pr_models=zeros(n_Patch,100);

pCount=0;
%select patch region and do model selection for each patch region
for(cPatch=1:49)
     
         testingSet=squeeze(testPatch(cPatch,:,:));
        
        gaussEst=patchLibrary{cPatch};
        
         
        %select testing set
        for(cPerson=1:size(imageIDTest,2))
            indices=find(imageIDTest(:,cPerson));
            l=length(indices);
            gallery(:,cPerson,cPatch)=testingSet(indices(1),:);
            probe(:,cPerson,cPatch)=testingSet(indices(2),:);
        end
end



        
;