function [allCaseD , allSTATS] = subXFStabSEG_v1_HY(BrainArea)


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

allCaseD = cell(3,3);

for gi = 1:3
    
    numCASES = sum(subTab.groupN == gi);
    
    caseINDS = find(subTab.groupN == gi);
    
    % Thalamus Volume
    volDif = nan(numCASES,1);
    
    % Normalized Volume
    normDif = nan(numCASES,1);
    
    % Up or Down in size
    signDir = nan(numCASES,1);
    
    
    ti = 1;
    
    for si = 1:numCASES
        
        cI = caseINDS(si);
        fsurgC = subTab.f_surg_n(cI);
        ssurgC = subTab.s_surg_n(cI);
        
        fsurgInd = ismember(totalSegTable.CaseName,fsurgC);
        ssurgInd = ismember(totalSegTable.CaseName,ssurgC);
        
        fsurgTab = totalSegTable(fsurgInd,:);
        ssurgTab = totalSegTable(ssurgInd,:);
        
        sideIND = subTab.f_surg_s{cI};
        
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
    
    allCaseD{1,gi} = volDif;
    allCaseD{2,gi} = normDif;
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


