function [allCaseD , allSTATS] = subXFStabHemi_v4_HY(ICVFlag)




bNames = {'entorhinal','lateralorbitofrontal','medialorbitofrontal'};

cd('C:\Users\johna\Dropbox\Publications_Meta\InProgress\Yilma_OB_DBS\DATA')

subTab = readtable('ob_subj_data.csv');

[allCaseD.PD] = getDATA(2, subTab, bNames, ICVFlag , 'PD');
[allSTATS.PD] = getSTATS(allCaseD.PD, 2, bNames);



end




function [dataOUTmain] = getDATA(groupNum, subTab, brFlag,  ICVFlag , condI)

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
        areaDif = nan(numCASES,1);
        
        ti = 1;
        
        for si = 1:numCASES
            
            cI = caseINDS(si);
            fsurgC = subTabt.f_surg_n(cI);
            ssurgC = subTabt.s_surg_n(cI);
            
            sideIND = subTabt.f_surg_s{cI};
            
            load('TotalSLHTable.mat');
            load('TotalSRHTable.mat');
            tableNL = totalLHTable;
            tableNL.CaseName = cellfun(@(x) str2double(x), tableNL.CaseName);
            tableNR = totalRHTable;
            tableNR.CaseName = cellfun(@(x) str2double(x), tableNR.CaseName);
            
            fsurgL = tableNL(ismember(tableNL.CaseName,fsurgC),:);
            fsurgR = tableNR(ismember(tableNR.CaseName,fsurgC),:);
            ssurgL = tableNL(ismember(tableNL.CaseName,ssurgC),:);
            ssurgR = tableNR(ismember(tableNR.CaseName,ssurgC),:);
            
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
            
            
            if isempty(fsurgL) || isempty(fsurgR) || isempty(ssurgL) || isempty(ssurgR)
                continue
            else
                
                if strcmp(sideIND,'L')
                    
                    FSindS = fsurgL.SurfArea{ismember(fsurgL.StructName,tmpBr)}/FSICV;
                    SSindS = ssurgL.SurfArea{ismember(ssurgL.StructName,tmpBr)}/SSICV;
                    
                    FSindN = fsurgR.SurfArea{ismember(fsurgR.StructName,tmpBr)}/FSICV;
                    SSindN = ssurgR.SurfArea{ismember(ssurgR.StructName,tmpBr)}/SSICV;
                    
                else
                    
                    FSindS = fsurgR.SurfArea{ismember(fsurgR.StructName,tmpBr)}/FSICV;
                    SSindS = ssurgR.SurfArea{ismember(ssurgR.StructName,tmpBr)}/SSICV;
                    
                    FSindN = fsurgL.SurfArea{ismember(fsurgL.StructName,tmpBr)}/FSICV;
                    SSindN = ssurgL.SurfArea{ismember(ssurgL.StructName,tmpBr)}/SSICV;
                    
                end
                
                FSdiff = (FSindS - FSindN) / max([FSindS , FSindN]);
                SSdiff = (SSindS - SSindN) / max([SSindS , SSindN]);
                
                areaDif(ti,1) = (FSdiff - SSdiff)*100;
                
                ti = ti + 1;
                
            end
        end
        
        areaDif = areaDif(~isnan(areaDif));
        
        dataOUTsub{1,gi} = areaDif;
        
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

