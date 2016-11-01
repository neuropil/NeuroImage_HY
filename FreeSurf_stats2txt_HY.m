function [] = FreeSurf_stats2txt_HY()



allFsLoc = 'Z:\Yilma_Project\Case_Data';
cd(allFsLoc);

dirFolds = dir;
dirFolds2 = {dirFolds.name};
dirFoldsA = dirFolds2(4:end);


for di = 1:length(dirFoldsA)
   
    
    tmpLoc = char(string(allFsLoc) + '\' + string(dirFoldsA{di}) + '\Freesurfer\stats');
    
    cd(tmpLoc)
    
    saveloc = 'Z:\Yilma_Project\CompiledCSVdata\'; 
    
    copyfile('aseg.stats',     fullfile(saveloc, [dirFoldsA{di} , '_SC.txt']),'f')
    copyfile('wmparc.stats',   fullfile(saveloc, [dirFoldsA{di} , '_WM.txt']),'f')
    copyfile('lh.aparc.stats', fullfile(saveloc, [dirFoldsA{di} , '_LH.txt']),'f')
    copyfile('rh.aparc.stats', fullfile(saveloc, [dirFoldsA{di} , '_RH.txt']),'f')
    
end








% allFsLoc = 'Z:\BRAiN_Project\All_Raw_Freesurfer';
% oldCodeLoc = 'E:\Dropbox\ForJeanelleCSV';