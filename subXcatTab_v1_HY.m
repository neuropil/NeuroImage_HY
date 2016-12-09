function [allCaseD , allSTATS] = subXcatTab_v1_HY(param)


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

allCaseD = cell(3,3);

for gi = 1:3
    
    numCASES = sum(subTab.groupN == gi);
    
    caseINDS = find(subTab.groupN == gi);
    
    % TIV
    absDif = nan(numCASES,1);
    
    % WM_MNI
    signDir = nan(numCASES,1);
    
    % WM_NAT
    relDif = nan(numCASES,1);
    
    ti = 1;
    
    for si = 1:numCASES
        
        cI = caseINDS(si);
        fsurgC = subTab.f_surg_n(cI);
        ssurgC = subTab.s_surg_n(cI);
        
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
    
    allCaseD{1,gi} = relDif;
    allCaseD{2,gi} = absDif;
    allCaseD{3,gi} = signDir;
    
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


