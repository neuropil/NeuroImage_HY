function [allCaseD , allSTATS] = subXFStabCortM_v3_HY()


cd('Z:\Yilma_Project\CompiledCSVdata')

subTab = readtable('ob_subj_data.csv');


            [allCaseD.PD] = getDATA(2, subTab, 'PD');
            [allSTATS.PD] = getSTATS(allCaseD.PD, 2);


end



function [dataOUT] = getDATA(groupNum, subTab, condI)

condIn = ismember(subTab.cond,condI);
subTabt = subTab(condIn,:);

dataOUT = cell(1,groupNum);

for gi = 1:groupNum
    
    numCASES = sum(subTabt.groupN == gi);
    
    caseINDS = find(subTabt.groupN == gi);
    
    % Thalamus Volume
    totalCM = nan(numCASES,1);
    
    
    ti = 1;
    
    for si = 1:numCASES
        
        cI = caseINDS(si);
        fsurgC = subTabt.f_surg_n(cI);
        ssurgC = subTabt.s_surg_n(cI);
        
        sideIND = subTabt.f_surg_s{cI};
        
        load('TotalBWMTable.mat');
        
        tableN = totalBrainWMTable;
        
        tableN.Case = cellfun(@(x) str2double(x), tableN.Case);
        
        fsurgInd = ismember(tableN.Case,fsurgC);
        ssurgInd = ismember(tableN.Case,ssurgC);
        
        fsurgTab = tableN(fsurgInd,:);
        ssurgTab = tableN(ssurgInd,:);
        
        if sum(fsurgInd) == 0 || sum(ssurgInd) == 0
            continue
        else
            
            if strcmp(sideIND,'L')
                brainFlagS = {'lhCorticalWhiteMatter,','lhCerebralWhiteMatter,'};
                brainFlagN = {'rhCorticalWhiteMatter,','rhCerebralWhiteMatter,'};
            else
                brainFlagS = {'rhCorticalWhiteMatter,','rhCerebralWhiteMatter,'};
                brainFlagN = {'lhCorticalWhiteMatter,','lhCerebralWhiteMatter,'};
            end
            
            
            
            FSIndS = fsurgTab.Volmm3(ismember(fsurgTab.BArea,brainFlagS));
            SSIndS = ssurgTab.Volmm3(ismember(ssurgTab.BArea,brainFlagS));
            
            FSIndN = fsurgTab.Volmm3(ismember(fsurgTab.BArea,brainFlagN));
            SSIndN = ssurgTab.Volmm3(ismember(ssurgTab.BArea,brainFlagN));
            
            
            FsDiff = (FSIndS - FSIndN) / max([FSIndN , FSIndS]);
            SsDiff = (SSIndS - SSIndN) / max([SSIndS , SSIndN]);
            
            perDif = (SsDiff - FsDiff)*100;
            
            
            
            
            
            totalCM(ti,1) = perDif;
            
            
            
            ti = ti + 1;
        end
    end
    
    totalCM = totalCM(~isnan(totalCM));
    
    
    dataOUT{1,gi} = totalCM;
    
    
end


end



function [statsOUT] = getSTATS(allCaseDpd, groupNum)

statsOUT = cell(1,groupNum);



forSTATS.data = [];
forSTATS.group = [];

for g2 = 1:groupNum
    
    forSTATS.data = [forSTATS.data ; allCaseDpd{1,g2}];
    
    grTmp = repmat(g2 , length(allCaseDpd{1,g2}) , 1);
    
    forSTATS.group = [forSTATS.group ; grTmp];
    
end

statsOUT{1,1} = forSTATS;




end
