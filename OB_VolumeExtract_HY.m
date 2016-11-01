function [] = OB_VolumeExtract_HY()


% First loop

fmainLoc = 'Z:\Yilma_Project\Case_Data';

cd(fmainLoc)

fdir1 = dir;
fdir2 = {fdir1.name};
fdir3 = fdir2(4:end);

fdir4 = fdir3(cellfun(@(x) ~isempty(strfind(x,'Case')), fdir3, 'UniformOutput',true));

for i = 1:length(fdir4)
    
    tmpD = [fmainLoc,'\',fdir4{i},'\NIFTI'];
    cd(tmpD);
    
    nameParts = strsplit(fdir4{i},'_');

    newName = ['ob_volumes_c',nameParts{2},'.csv'];
    
    newNloc = ['Z:\Yilma_Project\CompiledCSVdata\OV_CSV\',newName];   
    
    if ~exist([tmpD , '\ob_volumes.csv'],'file')
        continue
    else
        copyfile([tmpD , '\ob_volumes.csv'],newNloc);
    end
end



















end