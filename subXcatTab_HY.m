function [allCaseD , allSTATS] = subXcatTab_HY()


cd('Z:\Yilma_Project\CompiledCSVdata')

subTab = readtable('ob_subj_data.csv');

load('CAT12_All.mat');

allTable.CaseNum2 = cellfun(@(x) str2double(x), allTable.CaseNum); %#ok<NODEF>

allCaseD = cell(2,3);

for gi = 1:3
    
    numCASES = sum(subTab.groupN == gi);
    
    caseINDS = find(subTab.groupN == gi);
    
    % TIV
    diffTIV = nan(numCASES,1);
    
    % WM
    diffWMv = nan(numCASES,1);
    
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
            
            % TIV
            diffTIV(ti,1) = abs(double(allTable.TIV(ssurgInd) - allTable.TIV(fsurgInd))) / double(allTable.TIV(fsurgInd));
%             diffWMv(ti,1) = double(allTable.wm_nat(ssurgInd) - allTable.wm_nat(fsurgInd));
%             diffWMv(ti,1) = abs(double(allTable.wm_nat(ssurgInd) - allTable.wm_nat(fsurgInd))) / double(allTable.wm_nat(fsurgInd));
            
            diffWMv(ti,1) = abs(double(allTable.wm_mni(ssurgInd) - allTable.wm_mni(fsurgInd))) / double(allTable.wm_mni(fsurgInd));



            ti = ti + 1;
        end
    end
    
    diffWMv = diffWMv(~isnan(diffWMv));
    diffTIV = diffTIV(~isnan(diffTIV));
    allCaseD{1,gi} = diffTIV;
    allCaseD{2,gi} = diffWMv;
    
end


allSTATS = cell(1,2);

for si = 1:2
    
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


