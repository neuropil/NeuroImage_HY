function [missedFiles] = FreeSurf_hippoMove_HY(stage)


switch stage
    
    case 1
        
        allFsLoc = 'Z:\Yilma_Project\Case_Data';
        cd(allFsLoc);
        
        dirFolds = dir;
        dirFolds2 = {dirFolds.name};
        dirFoldsA = dirFolds2(4:end);
        
        missedFiles = cell(length(dirFoldsA),2);
        for di = 1:length(dirFoldsA)
            
            
            tmpLoc = char(string(allFsLoc) + '\' + string(dirFoldsA{di}) + '\Freesurfer\mri');
            
            if ~exist(tmpLoc , 'dir')
                continue
            end
            
            cd(tmpLoc)
            
            saveloc = 'Z:\Yilma_Project\CompiledCSVdata\Hippo_TXT';
            
            if ~exist('lh.hippoSfVolumes-T1.v10.txt' , 'file')
                missedFiles{di,1} = [dirFoldsA{di} , '_LHippo.txt'];
            else
                copyfile('lh.hippoSfVolumes-T1.v10.txt',   fullfile(saveloc, [dirFoldsA{di} , '_LHippo.txt']),'f')
            end
            
            if ~exist('rh.hippoSfVolumes-T1.v10.txt' , 'file')
                missedFiles{di,2} = [dirFoldsA{di} , '_RHippo.txt'];
            else
                copyfile('rh.hippoSfVolumes-T1.v10.txt',   fullfile(saveloc, [dirFoldsA{di} , '_RHippo.txt']),'f')
            end
            
            
            
        end
        
    case 2
        
        % Create master mat file
        
        txtLoc = 'Z:\Yilma_Project\CompiledCSVdata\Hippo_TXT';
        cd(txtLoc)
        
        dirList = dir('*.txt');
        
        flist = {dirList.name};
        
        iLines = cell(length(flist),1);
        
        hippoSegTable = table;
        
        for fi = 1:length(flist)
            
            nameParts = strsplit(flist{fi},{'_','.'});
            
            caseID = nameParts{2};
            sideID = nameParts{3};
            
            fid = fopen(flist{fi});
            tline = fgets(fid);
            allLines = cell(10000,1);
            lcount = 1;
            while ischar(tline)
                allLines{lcount,1} = tline;
                lcount = lcount + 1;
                tline = fgets(fid);
            end
            
            allLind = cellfun(@(x) ~isempty(x), allLines);
            allLines = allLines(allLind);
            iLines{fi} = allLines;

            [thippoTab] = extractFiles(allLines,caseID,sideID);
            hippoSegTable = [hippoSegTable ; thippoTab]; %#ok<*AGROW>

            
            
            
        end
        
        fclose(fid);
        missedFiles = nan;
        
end



cd('Z:\Yilma_Project\CompiledCSVdata')


% Subcortical Cortical
save('TotalHippoTable.mat','hippoSegTable');

writetable(hippoSegTable,'tHippoTable.csv')









end



function [brainTable] = extractFiles(inLINES, caseID, sideID)
%%%% cubic millimeters

dataNames = cell(length(inLINES),1);
dataVals = cell(length(inLINES),1);
caseCol = cell(length(inLINES),1);
sideCol = cell(length(inLINES),1);
colNames = {'Case','HippoArea','Volmm3','HemiS'};

for li = 1:length(inLINES)
    
    tparts = strsplit(inLINES{li});

    dataNames{li} = tparts{1};
    dataVals{li} = round(str2double(tparts{2}),3);
    caseCol{li} = caseID;
    
    if strcmp(sideID,'LHippo')
        sideCol{li} = 'L';
    else
        sideCol{li} = 'R';
    end

end

brainAind = cellfun(@(x) ~isempty(x), dataNames);
dataNames = dataNames(brainAind);
dataVals = dataVals(brainAind);

allCells = [caseCol, dataNames, dataVals , sideCol];
brainTable = cell2table(allCells,'VariableNames',colNames);

end

