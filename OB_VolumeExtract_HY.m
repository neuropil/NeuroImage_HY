function [] = OB_VolumeExtract_HY()


% First loop

fmainLoc = 'Z:\Yilma_Project\Case_Data';

cd(fmainLoc)

fdir1 = dir;
fdir2 = {fdir1.name};
fdir3 = fdir2(4:end);

fdir4 = fdir3(cellfun(@(x) ~isempty(strfind(x,'Case')), fdir3, 'UniformOutput',true));

for i = 1:length(fdir4)
    
    tmpD = [fmainLoc,'\',fdir4{i},'\NIFTI\OB_Volume'];
    cd(tmpD);
    
    fdir = dir('*.csv');
    fdirNs = {fdir.name};
    
    if isempty(fdirNs)
        continue
    end
    
    nameParts = strsplit(fdir4{i},'_');
    
    newName = ['ob_volumes_c',nameParts{2},'.csv'];
    
    newNloc = ['Z:\Yilma_Project\CompiledCSVdata\OV_CSV\',newName];
    
    if ~exist([tmpD , '\ob_volumes.csv'],'file') && ~exist([tmpD , '\ob_volumesC.csv'],'file')
        continue
    else
        if strcmp(fdirNs,'ob_volumes.csv')
            copyfile([tmpD , '\ob_volumes.csv'],newNloc);
        else
            copyfile([tmpD , '\ob_volumesC.csv'],newNloc);
        end
    end
end



















end