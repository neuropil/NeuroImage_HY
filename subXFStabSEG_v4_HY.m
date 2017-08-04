function [allCaseD , allSTATS] = subXFStabSEG_v4_HY(ICVFlag)


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

[allCaseD.PD] = getDATA(2, subTab, totalSegTable, bNames, ICVFlag ,'PD');
[allSTATS.PD] = getSTATS(allCaseD.PD, 2, bNames);

end






function [dataOUTmain] = getDATA(groupNum, subTab, totalSegTable, brFlag, ICVFlag, condI)

condIn = ismember(subTab.cond,condI);
subTabt = subTab(condIn,:);


load('CAT12_All.mat');
load('TotalBTable.mat');

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
            
            FS_ICVVal_CAT = double(allTable.TIV(ismember(allTable.CaseNum,num2str(fsurgC))));
            SS_ICVVal_CAT = double(allTable.TIV(ismember(allTable.CaseNum,num2str(ssurgC))));
            
            % Freesurfer
            FS_Tab_FS = totalBrainTable(ismember(totalBrainTable.Case,num2str(fsurgC)),:); %#ok<NODEF>
            SS_Tab_FS = totalBrainTable(ismember(totalBrainTable.Case,num2str(ssurgC)),:);
            
            FS_ICVVal_FS = double(FS_Tab_FS.Volmm3(ismember(FS_Tab_FS.BArea,'EstimatedTotalIntraCranialVol,')));
            SS_ICVVal_FS = double(SS_Tab_FS.Volmm3(ismember(SS_Tab_FS.BArea,'EstimatedTotalIntraCranialVol,')));
            
            if ICVFlag == 1
                
                FSICV = FS_ICVVal_CAT;
                SSICV = SS_ICVVal_CAT;
                
            elseif ICVFlag == 2
                
                FSICV = FS_ICVVal_FS;
                SSICV = SS_ICVVal_FS;
            end
            
            
            
            if isempty(fsurgTab) || isempty(ssurgTab)
                continue
            else
                
                if strcmp(sideIND,'L')
                    FSindS = fsurgTab.Volume_mm3{ismember(fsurgTab.StructName,['Left',tmpBr])}/FSICV;
                    SSindS = ssurgTab.Volume_mm3{ismember(ssurgTab.StructName,['Left',tmpBr])}/SSICV;
                    
                    FSindN = fsurgTab.Volume_mm3{ismember(fsurgTab.StructName,['Right',tmpBr])}/FSICV;
                    SSindN = ssurgTab.Volume_mm3{ismember(ssurgTab.StructName,['Right',tmpBr])}/SSICV;
                    
                else
                    
                    FSindS = fsurgTab.Volume_mm3{ismember(fsurgTab.StructName,['Right',tmpBr])}/FSICV;
                    SSindS = ssurgTab.Volume_mm3{ismember(ssurgTab.StructName,['Right',tmpBr])}/SSICV;
                    
                    FSindN = fsurgTab.Volume_mm3{ismember(fsurgTab.StructName,['Left',tmpBr])}/FSICV;
                    SSindN = ssurgTab.Volume_mm3{ismember(ssurgTab.StructName,['Left',tmpBr])}/SSICV;
                    
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

