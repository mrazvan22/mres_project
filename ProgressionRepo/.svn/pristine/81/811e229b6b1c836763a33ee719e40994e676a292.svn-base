function [destinationDir,patchLibrary, testingSet, imageIDTest, name]=dataSetsTest(dataType,N_T_DIST,modelType)

switch dataType
    
    case 1 % ====== data with ground truth illumination XM2VTS =========
            if (modelType==1)
            libDir='C:\PhD Research\Face Recognition\MGM\Face Mix Gauss\mixT\Patch Version\Standard Model\Mat Data\XM2VTS\';
            libName=['XM2VTSLibraryVersion2T' num2str(N_T_DIST) '.mat'];
            destinationDir='C:\PhD Research\Face Recognition\MGM\Face Mix Gauss\mixT\Patch Version\Standard Model\Mat Data\XM2VTS\';
            name=['XM2VTSVersion2T' num2str(N_T_DIST) '.mat'];
            
            else
            libDir='C:\PhD Research\Face Recognition\MGM\Face Mix Gauss\mixT\Patch Version\Appearance Model\Mat Data\XM2VTS\';
            libName=['XM2VTSLibraryVersion3T' num2str(N_T_DIST) '.mat'];
            destinationDir='C:\PhD Research\Face Recognition\MGM\Face Mix Gauss\mixT\Patch Version\Appearance Model\Mat Data\XM2VTS\';
            name=['XM2VTSVersion3T' num2str(N_T_DIST) '.mat'];
            end
            
            %load the library
            load([libDir libName]);
            
            dataDir='C:\PhD Research\Face Recognition\MGM\Face Mix Gauss\Version 2\Mat data\';
            load ([dataDir 'patchTestingSets8.mat'])
            testingSet=testPatch;
            %Load the dataIDmatrix
            sourceDir='C:\PhD Research\Face Recognition\Probabilistic FR\Enhanced Images\'
            load ([sourceDir 'enhancedSets.mat'],'imageIDTest');
            
            

        
    case 2 % Gabor filtered data with ground truth illumination XM2VTS
        
    case 3 % light data XM2VTS
    
    case 4 % Gabor filtered light data XM2VTS 
end