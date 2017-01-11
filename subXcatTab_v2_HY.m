function [allCaseD , allSTATS] = subXcatTab_v2_HY(param)


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
condS = {'PD','ET'};
for ci = 1:2
    switch ci
        case 1
            [allCaseD.PD] = getDATA(3, subTab, allTable, pUSE, condS{ci});
            [allSTATS.PD] = getSTATS(allCaseD.PD, 3);
        case 2
            [allCaseD.ET] = getDATA(1, subTab, allTable, pUSE, condS{ci});
            [allSTATS.ET] = getSTATS(allCaseD.ET, 1);
    end
end



end



function [dataOUT] = getDATA(groupNum, subTab, allTable, pUSE, condI)

condIn = ismember(subTab.cond,condI);
subTabt = subTab(condIn,:);

dataOUT = cell(3,groupNum);

for gi = 1:groupNum
    
    numCASES = sum(subTabt.groupN == gi);
    
    caseINDS = find(subTabt.groupN == gi);
    
    % TIV
    absDif = nan(numCASES,1);
    
    % WM_MNI
    signDir = nan(numCASES,1);
    
    % WM_NAT
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

            relDif(ti,1) = abs(double(allTable.(pUSE)(ssurgInd) - allTable.(pUSE)(fsurgInd))) / double(allTable.(pUSE)(fsurgInd));

            absDif(ti,1) = double(allTable.(pUSE)(ssurgInd) - allTable.(pUSE)(fsurgInd)) / double(allTable.(pUSE)(fsurgInd));
            
            signDir(ti,1) = sign(absDif(ti,1));
            
            ti = ti + 1;
        end
    end
    
    relDif = relDif(~isnan(relDif));
    absDif = absDif(~isnan(absDif));
    signDir = signDir(~isnan(signDir));
    
    dataOUT{1,gi} = relDif;
    dataOUT{2,gi} = absDif;
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

