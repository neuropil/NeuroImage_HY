function [allCaseD , allSTATS] = subXFStabSEG_v3_HY()


% switch BrainArea
%     case 'Thal'
%
%         brFlag = '-Thalamus-Proper';
%
%     case 'Caud'
%
%         brFlag = '-Caudate';
%
%     case 'Hippo'
%
%         brFlag = '-Hippocampus';
% end

bNames = {'-Thalamus-Proper','-Caudate'};

cd('Z:\Yilma_Project\CompiledCSVdata')

subTab = readtable('ob_subj_data.csv');

load('TotalSTable.mat');

totalSegTable.CaseName = cellfun(@(x) str2double(x), totalSegTable.CaseName); %#ok<NODEF>

[allCaseD.PD] = getDATA(2, subTab, totalSegTable, bNames,'PD');
[allSTATS.PD] = getSTATS(allCaseD.PD, 2, bNames);

end






function [dataOUTmain] = getDATA(groupNum, subTab, totalSegTable, brFlag, condI)

condIn = ismember(subTab.cond,condI);
subTabt = subTab(condIn,:);

dataOUTmain = cell(length(brFlag),1);
for aai = 1:length(brFlag)
    
    tmpBr = brFlag{aai};
    
    dataOUTsub = cell(1,groupNum);
    
    for gi = 1:groupNum
        
        
        numCASES = sum(subTabt.groupN == gi);
        
        caseINDS = find(subTabt.groupN == gi);
        
        % Thalamus Volume
        volDif = nan(numCASES,1);
        
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
            if isempty(fsurgTab) || isempty(ssurgTab) 
                continue
            else
                
                if strcmp(sideIND,'L')
                    FSindS = fsurgTab.Volume_mm3{ismember(fsurgTab.StructName,['Left',tmpBr])};
                    SSindS = ssurgTab.Volume_mm3{ismember(ssurgTab.StructName,['Left',tmpBr])};
                    
                    FSindN = fsurgTab.Volume_mm3{ismember(fsurgTab.StructName,['Right',tmpBr])};
                    SSindN = ssurgTab.Volume_mm3{ismember(ssurgTab.StructName,['Right',tmpBr])};
                    
                else
                    
                    FSindS = fsurgTab.Volume_mm3{ismember(fsurgTab.StructName,['Right',tmpBr])};
                    SSindS = ssurgTab.Volume_mm3{ismember(ssurgTab.StructName,['Right',tmpBr])};
                    
                    FSindN = fsurgTab.Volume_mm3{ismember(fsurgTab.StructName,['Left',tmpBr])};
                    SSindN = ssurgTab.Volume_mm3{ismember(ssurgTab.StructName,['Left',tmpBr])};
                    
                end
                
                
                
                
                FSdiff = (FSindS - FSindN) / max([FSindS , FSindN]);
                SSdiff = (SSindS - SSindN) / max([SSindS , SSindN]);
                
                volDif(ti,1) = (FSdiff - SSdiff)*100;
                
                ti = ti + 1;
            end
        end
        
        volDif = volDif(~isnan(volDif));
        
        
        dataOUTsub{1,gi} = volDif;
    end
    dataOUTmain{aai} = dataOUTsub;
    
end


end









function [statsOUT] = getSTATS(allCaseDpd, groupNum, brNames)

statsOUT = cell(length(allCaseDpd),1);

for si = 1:length(allCaseDpd)
    
    forSTATS.data = [];
    forSTATS.group = [];
    forSTATS.brN = brNames{si};
    
    for g2 = 1:groupNum
        
        forSTATS.data = [forSTATS.data ; allCaseDpd{si}{g2}];
        
        grTmp = repmat(g2 , length(allCaseDpd{si}{g2}) , 1);
        
        forSTATS.group = [forSTATS.group ; grTmp];
        
    end
    
    statsOUT{si,1} = forSTATS;
    
end



end

