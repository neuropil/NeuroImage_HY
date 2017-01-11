function [allCaseD , allSTATS] = subXFStabHemi_v2_HY(BrainArea)


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



condS = {'PD','ET'};
for ci = 1:2
    switch ci
        case 1
            [allCaseD.PD] = getDATA(3, subTab, brFlag, condS{ci});
            [allSTATS.PD] = getSTATS(allCaseD.PD, 3);
        case 2
            [allCaseD.ET] = getDATA(1, subTab, brFlag, condS{ci});
            [allSTATS.ET] = getSTATS(allCaseD.ET, 1);
    end
end



end




function [dataOUT] = getDATA(groupNum, subTab, brFlag, condI)

condIn = ismember(subTab.cond,condI);
subTabt = subTab(condIn,:);

dataOUT = cell(3,groupNum);

for gi = 1:groupNum
    
    numCASES = sum(subTabt.groupN == gi);
    
    caseINDS = find(subTabt.groupN == gi);
    
    % Thalamus Volume
    areaDif = nan(numCASES,1);
    
    % Normalized Volume
    volDif = nan(numCASES,1);
    
    % Up or Down in size
    thickAve = nan(numCASES,1);
    
    
    ti = 1;
    
    for si = 1:numCASES
        
        cI = caseINDS(si);
        fsurgC = subTabt.f_surg_n(cI);
        ssurgC = subTabt.s_surg_n(cI);
        
        sideIND = subTabt.f_surg_s{cI};
        
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
    
    dataOUT{1,gi} = areaDif;
    dataOUT{2,gi} = volDif;
    dataOUT{3,gi} = thickAve;
    
end


end











function [statsOUT] = getSTATS(allCaseDpd, groupNum)

statsOUT = cell(1,groupNum);

for si = 1:3
    
    forSTATS.data = [];
    forSTATS.group = [];
    
    for g2 = 1:groupNum
        
        forSTATS.data = [forSTATS.data ; allCaseDpd{si,g2}];
        
        grTmp = repmat(g2 , length(allCaseDpd{si,g2}) , 1);
        
        forSTATS.group = [forSTATS.group ; grTmp];
        
    end
    
    statsOUT{1,si} = forSTATS;
    
end



end

