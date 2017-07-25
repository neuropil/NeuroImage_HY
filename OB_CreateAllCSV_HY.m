function [] = OB_CreateAllCSV_HY()

% Location

mainVolLoc = 'Z:\Yilma_Project\CompiledCSVdata\OV_CSV';

cd(mainVolLoc)

csvDir = dir('*.csv');
csvDirA = {csvDir.name};

allOBTab = table;
for ci = 1:length(csvDirA)
    
    % tmpTable
    tmpTable = readtable(csvDirA{ci});
    % tmpCname
    nParts = strsplit(csvDirA{ci},{'_','.'});
    tcName = nParts{3};
    tmpTable = tmpTable(tmpTable.LabelId ~= 0,:);
    
    % Label List

    labelID = {'R_OB';'L_OB'};
    
    hemiID = {'R';'L'};
    
    labelColor = [[255 0 0] ; [0 255 0]]/255;
    
    fslColors = cell(size(labelColor,1),1);
    for li = 1:size(labelColor,1)
        fslColors{li,1} = labelColor(li,:);
    end

        tmpTable.LabelName = labelID;
        tmpTable.HemiN = hemiID;
    
    tmpTable.LabColor = fslColors;
    
    tmpTable.Properties.VariableNames([3 4 5 6]) = {'NVoxels' , 'Volume_mm3' , 'IMmean' , 'IMstddev'};
    
    PDcases = {'c205','c215','c217','c219','c220','c221','c222','c223','c224','c226','c227','c228','c229','c233',...
        'c234','c235','c250','c251','c252','c253','c254','c255','c256','c258','c259','c262','c263','c264','c265','c266',...
        'c267','c268','c271','c272','c273','c274','c278','c283','c284','c287','c288','c289','c290','c291','c292','c293','c294',...
        'c295','c296','c297','c299','c300','c301','c303','c305','c312','c316','c326','c328','c330','c337','c341','c342','c343','c345','c358'};
    
    ETcases = {'c241','c249','c248','c270','c277','c286'};
    
    if ismember(tcName,PDcases)
        tmpTable.Condition = cellstr(repmat('PD',height(tmpTable),1));
    elseif ismember(tcName,ETcases)
        tmpTable.Condition = cellstr(repmat('ET',height(tmpTable),1));
    else
        keyboard
    end
    
    caseCol = cellstr(repmat(tcName,height(tmpTable),1));
    
    tmpTable.caseID = caseCol;
    
    
    
    
    
    
    allOBTab = [allOBTab ; tmpTable]; %#ok<AGROW>
    
end

cd('Z:\Yilma_Project\CompiledCSVdata');
% Subcortical Cortical
save('AllOBdata.mat','allOBTab');
writetable(allOBTab,'AllOBdata.csv')




end
