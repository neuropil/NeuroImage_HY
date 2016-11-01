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
    
    PDcases = {'c250','c251','c253','c258','c259','c263','c264','c265','c266',...
        'c271','c272','c274','c283','c284','c288','c289','c291',...
        'c296'};
    
    ETcases = {'c241','c249','c270','c277','c286'};
    
    if ismember(tcName,PDcases)
        tmpTable.Condition = cellstr(repmat('PD',height(tmpTable),1));
    elseif ismember(tcName,ETcases)
        tmpTable.Condition = cellstr(repmat('ET',height(tmpTable),1));
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
