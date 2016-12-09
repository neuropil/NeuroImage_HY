function [allCaseD , allSTATS] = subOBVOL_v1_HY(absFlag)

cd('Z:\Yilma_Project\CompiledCSVdata')

subTab = readtable('ob_subj_data.csv');

load('AllOBdata.mat');

allOBTab.caseID = cellfun(@(x) str2double(x(2:4)), allOBTab.caseID); %#ok<NODEF>

allCaseD = cell(2,3);

for gi = 1:3
    
    numCASES = sum(subTab.groupN == gi);
    
    caseINDS = find(subTab.groupN == gi);
    
    % Thalamus Volume
    volDif = nan(numCASES,1);
    
    % Normalized Volume
    normDif = nan(numCASES,1);
    
    ti = 1;
    
    for si = 1:numCASES
        
        cI = caseINDS(si);
        fsurgC = subTab.f_surg_n(cI);
        ssurgC = subTab.s_surg_n(cI);
        
        fsurgInd = ismember(allOBTab.caseID,fsurgC);
        ssurgInd = ismember(allOBTab.caseID,ssurgC);
        
        fsurgTab = allOBTab(fsurgInd,:);
        ssurgTab = allOBTab(ssurgInd,:);
        
        sideIND = subTab.f_surg_s{cI};
        
        if strcmp(sideIND,'L')
            thalFlag = 'L_OB';
        else
            thalFlag = 'R_OB';
        end
        
        if sum(fsurgInd) == 0 || sum(ssurgInd) == 0
            continue
        else
            
            FSindT = ismember(fsurgTab.LabelName,thalFlag);
            SSindT = ismember(ssurgTab.LabelName,thalFlag);
            
            if absFlag
                volDif(ti,1) = abs(double(ssurgTab.Volume_mm3(SSindT) - fsurgTab.Volume_mm3(FSindT))) / double(fsurgTab.Volume_mm3(FSindT));
            else
                volDif(ti,1) = double(ssurgTab.Volume_mm3(SSindT) - fsurgTab.Volume_mm3(FSindT)) / double(fsurgTab.Volume_mm3(FSindT));
            end
            normDif(ti,1) = abs(double(ssurgTab.IMmean(SSindT) - fsurgTab.IMmean(FSindT))) / double(fsurgTab.IMmean(FSindT));
            
            ti = ti + 1;
        end
    end
    
    volDif = volDif(~isnan(volDif));
    normDif = normDif(~isnan(normDif));
    
    
    allCaseD{1,gi} = volDif;
    allCaseD{2,gi} = normDif;
    
    
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


