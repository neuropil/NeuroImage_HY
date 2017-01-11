function [allCaseD , allSTATS] = subXFStabSEG_v2_HY(BrainArea)


switch BrainArea
    case 'Thal'
        
        brFlag = '-Thalamus-Proper';
        
    case 'Caud'
        
        brFlag = '-Caudate';
        
    case 'Hippo'
        
        brFlag = '-Hippocampus';
end


cd('Z:\Yilma_Project\CompiledCSVdata')

subTab = readtable('ob_subj_data.csv');

load('TotalSTable.mat');

totalSegTable.CaseName = cellfun(@(x) str2double(x), totalSegTable.CaseName); %#ok<NODEF>



condS = {'PD','ET'};
for ci = 1:2
    switch ci
        case 1
            [allCaseD.PD] = getDATA(3, subTab, totalSegTable, brFlag, condS{ci});
            [allSTATS.PD] = getSTATS(allCaseD.PD, 3);
        case 2
            [allCaseD.ET] = getDATA(1, subTab, totalSegTable, brFlag, condS{ci});
            [allSTATS.ET] = getSTATS(allCaseD.ET, 1);
    end
end




end






function [dataOUT] = getDATA(groupNum, subTab, totalSegTable, brFlag, condI)

condIn = ismember(subTab.cond,condI);
subTabt = subTab(condIn,:);

dataOUT = cell(3,groupNum);

for gi = 1:groupNum
    
    numCASES = sum(subTabt.groupN == gi);
    
    caseINDS = find(subTabt.groupN == gi);
    
    % Thalamus Volume
    volDif = nan(numCASES,1);
    
    % Normalized Volume
    normDif = nan(numCASES,1);
    
    % Up or Down in size
    signDir = nan(numCASES,1);
    
    
    ti = 1;
    
    for si = 1:numCASES
        
        cI = caseINDS(si);
        fsurgC = subTabt.f_surg_n(cI);
        ssurgC = subTabt.s_surg_n(cI);
        
        fsurgInd = ismember(totalSegTable.CaseName,fsurgC);
        ssurgInd = ismember(totalSegTable.CaseName,ssurgC);
        
        fsurgTab = totalSegTable(fsurgInd,:);
        ssurgTab = totalSegTable(ssurgInd,:);
        
        sideIND = subTabt.f_surg_s{cI};
        
        if strcmp(sideIND,'L')
            thalFlag = ['Left',brFlag];
        else
            thalFlag = ['Right',brFlag];
        end
        
        FSindT = ismember(fsurgTab.StructName,thalFlag);
        SSindT = ismember(ssurgTab.StructName,thalFlag);
        
        if sum(fsurgInd) == 0 || sum(ssurgInd) == 0
            continue
        else
            
            volDif(ti,1) = double(ssurgTab.Volume_mm3{SSindT} - fsurgTab.Volume_mm3{FSindT}) / double(fsurgTab.Volume_mm3{FSindT});
            
            normDif(ti,1) = abs(double(ssurgTab.normMean{SSindT} - fsurgTab.normMean{FSindT})) / double(fsurgTab.normMean{FSindT});
            
            signDir(ti,1) = sign(volDif(ti,1));
            
            ti = ti + 1;
        end
    end
    
    volDif = volDif(~isnan(volDif));
    normDif = normDif(~isnan(normDif));
    signDir = signDir(~isnan(signDir));
    
    dataOUT{1,gi} = volDif;
    dataOUT{2,gi} = normDif;
    dataOUT{3,gi} = signDir;
    
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

