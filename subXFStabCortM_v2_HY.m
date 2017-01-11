function [allCaseD , allSTATS] = subXFStabCortM_v2_HY()


cd('Z:\Yilma_Project\CompiledCSVdata')

subTab = readtable('ob_subj_data.csv');

condS = {'PD','ET'};
for ci = 1:2
    switch ci
        case 1
            [allCaseD.PD] = getDATA(3, subTab, condS{ci});
            [allSTATS.PD] = getSTATS(allCaseD.PD, 3);
        case 2
            [allCaseD.ET] = getDATA(1, subTab, condS{ci});
            [allSTATS.ET] = getSTATS(allCaseD.ET, 1);
    end
end


end



function [dataOUT] = getDATA(groupNum, subTab, condI)

condIn = ismember(subTab.cond,condI);
subTabt = subTab(condIn,:);

dataOUT = cell(3,groupNum);

for gi = 1:groupNum
    
    numCASES = sum(subTabt.groupN == gi);
    
    caseINDS = find(subTabt.groupN == gi);
    
    % Thalamus Volume
    totalCM = nan(numCASES,1);
    
    % Normalized Volume
    hemiCM = nan(numCASES,1);
    
    % Up or Down in size
    tiv = nan(numCASES,1);
    
    
    ti = 1;
    
    for si = 1:numCASES
        
        cI = caseINDS(si);
        fsurgC = subTabt.f_surg_n(cI);
        ssurgC = subTabt.s_surg_n(cI);
        
        sideIND = subTabt.f_surg_s{cI};
        
        if strcmp(sideIND,'L')
            brainFlag = 'lhCorticalWhiteMatter,';
        else
            brainFlag = 'rhCorticalWhiteMatter,';
        end
        
        load('TotalBWMTable.mat');

        tableN = totalBrainWMTable;
        
        tableN.Case = cellfun(@(x) str2double(x), tableN.Case); 
        
        fsurgInd = ismember(tableN.Case,fsurgC);
        ssurgInd = ismember(tableN.Case,ssurgC);
        
        fsurgTab = tableN(fsurgInd,:);
        ssurgTab = tableN(ssurgInd,:);

        FSindT1 = ismember(fsurgTab.BArea,brainFlag);
        SSindT1 = ismember(ssurgTab.BArea,brainFlag);
        
        FSindT2 = ismember(fsurgTab.BArea,'CorticalWhiteMatter,');
        SSindT2 = ismember(ssurgTab.BArea,'CorticalWhiteMatter,');
        
        FSindT3 = ismember(fsurgTab.BArea,'EstimatedTotalIntraCranialVol,');
        SSindT3 = ismember(ssurgTab.BArea,'EstimatedTotalIntraCranialVol,');
        
        if sum(fsurgInd) == 0 || sum(ssurgInd) == 0
            continue
        else
            
            totalCM(ti,1) = abs(double(ssurgTab.Volmm3(SSindT2) - fsurgTab.Volmm3(FSindT2))) / double(fsurgTab.Volmm3(FSindT2));
            
            hemiCM(ti,1) = abs(double(ssurgTab.Volmm3(SSindT1) - fsurgTab.Volmm3(FSindT1))) / double(fsurgTab.Volmm3(FSindT1));
            
            tiv(ti,1) = abs(double(ssurgTab.Volmm3(SSindT3) - fsurgTab.Volmm3(FSindT3))) / double(fsurgTab.Volmm3(FSindT3));
            
            ti = ti + 1;
        end
    end
    
    totalCM = totalCM(~isnan(totalCM));
    hemiCM = hemiCM(~isnan(hemiCM));
    tiv = tiv(~isnan(tiv));
    
    dataOUT{1,gi} = totalCM;
    dataOUT{2,gi} = hemiCM;
    dataOUT{3,gi} = tiv;
    
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
