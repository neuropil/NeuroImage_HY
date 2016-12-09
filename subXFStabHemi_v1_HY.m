function [allCaseD , allSTATS] = subXFStabHemi_v1_HY(BrainArea)


switch BrainArea
    case 'entor'
        
        brFlag = 'entorhinal';
        
    case 'latOB'
        
        brFlag = 'lateralorbitofrontal';
        
    case 'medOB'
        
        brFlag = 'medialorbitofrontal';
end


cd('Z:\Yilma_Project\CompiledCSVdata')

subTab = readtable('ob_subj_data.csv');



allCaseD = cell(3,3);

for gi = 1:3
    
    numCASES = sum(subTab.groupN == gi);
    
    caseINDS = find(subTab.groupN == gi);
    
    % Thalamus Volume
    areaDif = nan(numCASES,1);
    
    % Normalized Volume
    volDif = nan(numCASES,1);
    
    % Up or Down in size
    thickAve = nan(numCASES,1);
    
    
    ti = 1;
    
    for si = 1:numCASES
        
        cI = caseINDS(si);
        fsurgC = subTab.f_surg_n(cI);
        ssurgC = subTab.s_surg_n(cI);
        
        sideIND = subTab.f_surg_s{cI};
        
        if strcmp(sideIND,'L')
            toLoadfl = 'TotalSLHTable.mat';
             load(toLoadfl);
            tableN = totalLHTable;
        else
            toLoadfl = 'TotalSRHTable.mat';
             load(toLoadfl);
             tableN = totalRHTable;
        end

        tableN.CaseName = cellfun(@(x) str2double(x), tableN.CaseName); 
        
        fsurgInd = ismember(tableN.CaseName,fsurgC);
        ssurgInd = ismember(tableN.CaseName,ssurgC);
        
        fsurgTab = tableN(fsurgInd,:);
        ssurgTab = tableN(ssurgInd,:);

        FSindT = ismember(fsurgTab.StructName,brFlag);
        SSindT = ismember(ssurgTab.StructName,brFlag);
        
        if sum(fsurgInd) == 0 || sum(ssurgInd) == 0
            continue
        else
            
            areaDif(ti,1) = abs(double(ssurgTab.SurfArea{SSindT} - fsurgTab.SurfArea{FSindT})) / double(fsurgTab.SurfArea{FSindT});
            
            volDif(ti,1) = abs(double(ssurgTab.GrayVol{SSindT} - fsurgTab.GrayVol{FSindT})) / double(fsurgTab.GrayVol{FSindT});
            
            thickAve(ti,1) = abs(double(ssurgTab.ThickAvg{SSindT} - fsurgTab.ThickAvg{FSindT})) / double(fsurgTab.ThickAvg{FSindT});
            
            ti = ti + 1;
        end
    end
    
    areaDif = areaDif(~isnan(areaDif));
    volDif = volDif(~isnan(volDif));
    thickAve = thickAve(~isnan(thickAve));
    
    allCaseD{1,gi} = areaDif;
    allCaseD{2,gi} = volDif;
    allCaseD{3,gi} = thickAve;
    
end


allSTATS = cell(1,3);

for si = 1:3
    
    forSTATS.data = [];
    forSTATS.group = [];
    
    for g2 = 1:3
        
        forSTATS.data = [forSTATS.data ; allCaseD{si,g2}];
        
        grTmp = repmat(g2 , length(allCaseD{si,g2}) , 1);
        
        forSTATS.group = [forSTATS.group ; grTmp];
        
    end
    
    allSTATS{1,si} = forSTATS;
    
end




end


