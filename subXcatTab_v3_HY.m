function [allCaseD , allSTATS] = subXcatTab_v3_HY(param)


switch param
    case 'WMnat'
        pUSE = 'wm_nat';
        
    case 'WMmni'
        pUSE = 'wm_mni';
        
    case 'TIV'
        pUSE = 'TIV';
        
    case 'GMnat'
        pUSE = 'gm_nat';
        
    case 'GMmni'
        pUSE = 'gm_mni';
        
    case 'CTthick'
        pUSE = 'cortThick_ave';
end


cd('Z:\Yilma_Project\CompiledCSVdata')

subTab = readtable('ob_subj_data.csv');

load('CAT12_All.mat');

allTable.CaseNum2 = cellfun(@(x) str2double(x), allTable.CaseNum); %#ok<NODEF>

% allCaseDpd = cell(3,3);

[allCaseD.PD] = getDATA(2, subTab, allTable, pUSE, 'PD');
[allSTATS.PD] = getSTATS(allCaseD.PD, 2);



end



function [dataOUT] = getDATA(groupNum, subTab, allTable, pUSE, condI)

condIn = ismember(subTab.cond,condI);
subTabt = subTab(condIn,:);

dataOUT = cell(2,groupNum);

for gi = 1:groupNum
    
    numCASES = sum(subTabt.groupN == gi);
    
    caseINDS = find(subTabt.groupN == gi);
    
    % 
    absDif = nan(numCASES,1);
    
    % 
    relDif = nan(numCASES,1);
    
    ti = 1;
    
    for si = 1:numCASES
        
        cI = caseINDS(si);
        fsurgC = subTabt.f_surg_n(cI);
        ssurgC = subTabt.s_surg_n(cI);
        
        fsurgInd = ismember(allTable.CaseNum2,fsurgC);
        ssurgInd = ismember(allTable.CaseNum2,ssurgC);
        
        if sum(fsurgInd) == 0 || sum(ssurgInd) == 0
            continue
        else

%             relDif(ti,1) = (double(allTable.(pUSE)(fsurgInd)) - double(allTable.(pUSE)(ssurgInd))) / double(max([allTable.(pUSE)(fsurgInd) , allTable.(pUSE)(ssurgInd)]));
            relDif(ti,1) = double(allTable.(pUSE)(ssurgInd)) - double(allTable.(pUSE)(fsurgInd));
            absDif(ti,1) = double(allTable.(pUSE)(ssurgInd) - allTable.(pUSE)(fsurgInd)) / double(allTable.(pUSE)(fsurgInd));
            
            
            ti = ti + 1;
        end
    end
    
    relDif = relDif(~isnan(relDif));
    absDif = absDif(~isnan(absDif));

    
    dataOUT{1,gi} = relDif;
    dataOUT{2,gi} = absDif;

    
end


end



function [statsOUT] = getSTATS(allCaseDpd, groupNum)

statsOUT = cell(1,groupNum);

for si = 1:2
    
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

