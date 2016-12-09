function [allCaseD , allSTATS] = subXFStabCortM_v1_HY()


cd('Z:\Yilma_Project\CompiledCSVdata')

subTab = readtable('ob_subj_data.csv');



allCaseD = cell(3,3);

for gi = 1:3
    
    numCASES = sum(subTab.groupN == gi);
    
    caseINDS = find(subTab.groupN == gi);
    
    % Thalamus Volume
    totalCM = nan(numCASES,1);
    
    % Normalized Volume
    hemiCM = nan(numCASES,1);
    
    % Up or Down in size
    tiv = nan(numCASES,1);
    
    
    ti = 1;
    
    for si = 1:numCASES
        
        cI = caseINDS(si);
        fsurgC = subTab.f_surg_n(cI);
        ssurgC = subTab.s_surg_n(cI);
        
        sideIND = subTab.f_surg_s{cI};
        
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
    
    allCaseD{1,gi} = totalCM;
    allCaseD{2,gi} = hemiCM;
    allCaseD{3,gi} = tiv;
    
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


