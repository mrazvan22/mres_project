function rt=buildImage(patchNum1,patchNum2,reconsPatch,patch_sizeX,patch_sizeY,nZ)

a=[1:patchNum1*patchNum2];
a=reshape(a,patchNum1,patchNum2)';

for(i=1:patchNum1*patchNum2)
    i
    patchMatrix{i}=reshape(reconsPatch(:,i),patch_sizeX,patch_sizeY,nZ);
end
for(i=1:patchNum1)
img{i}=cat(1,patchMatrix{a(:,i)});
end

reconX=cat(2,img{:});

rt=reconX;