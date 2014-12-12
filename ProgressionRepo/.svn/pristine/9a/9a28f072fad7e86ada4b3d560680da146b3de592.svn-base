function rt=createConseqDataID(cSample, cData, cIndiv)

destinationDir='C:\PhD Research\Face Recognition\Probabilistic FR\Version3\';
column=zeros(cData,1);
column(1:cSample)=1;

count=0;
dataIDMatrix(:,1)=column;

for(i=2:cIndiv)
    dataIDMatrix(:,i)=circshift(column,i+count);
    count=count+1;
end

%save([destinationDir 'testDataIDMatrix.mat'],'dataIDMatrix');

rt=dataIDMatrix;
