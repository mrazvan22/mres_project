function [trainingSet, imageIDTrain, name]=dataSetsTrain(dataType)

switch dataType
    
    case 1 % ====== data with ground truth illumination XM2VTS =========
            libDir='C:\PhD Research\Face Recognition\MGM\Face Mix Gauss\Version 2\Mat data\';
            load ([libDir 'patchTrainingSets8.mat'])
            trainingSet=trainPatch;
            %Load the dataIDmatrix
            sourceDir='C:\PhD Research\Face Recognition\Probabilistic FR\Enhanced Images\'
            load ([sourceDir 'enhancedSets.mat'],'imageIDTrain');
            name='XM2VTSLibrary';
            

        
    case 2 % Gabor filtered data with ground truth illumination XM2VTS
        
    case 3 % light data XM2VTS
    
    case 4 % Gabor filtered light data XM2VTS 
end