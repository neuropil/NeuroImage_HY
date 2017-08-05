function [allCaseD , allSTATS ,groupIDS] = subXhippo_HYb(outFlag,ICVFlag,comp)


if comp == 1
    
    cd('Z:\Yilma_Project\CompiledCSVdata')
    
elseif comp == 2
    
    cd('C:\Users\johna\Dropbox\Publications_Meta\InProgress\Yilma_OB_DBS\DATA')
    
end

subTab = readtable('ob_subj_data.csv');

load('TotalHippoTable.mat');

hippoSegTable.Case = cellfun(@(x) str2double(x), hippoSegTable.Case); %#ok<NODEF>

% bNames = unique(hippoSegTable.HippoArea);

bNames = {'CA1','CA3','CA4','Whole_hippocampus'};

[allCaseD, groupIDS] = getDATA(2, subTab, hippoSegTable, bNames, ICVFlag,outFlag, 'PD');
[allSTATS.PD] = getSTATS(allCaseD, 2 , bNames);




end






function [dataOUTmain, groupIDS] = getDATA(groupNum, subTab, totalSegTable, brFlag, ICVFlag,outFlag, condI)

condIn = ismember(subTab.cond,condI);
subTabt = subTab(condIn,:);

dataOUTmain = cell(length(brFlag),1);



load('CAT12_All.mat');
load('TotalBTable.mat');

for aai = 1:length(brFlag)
    
    tmpBr = brFlag{aai};
    
    dataOUTsub = cell(2,groupNum);
    
    groupIDS = cell(2,groupNum);
    
    for gi = 1:groupNum
        
        
        numCASES = sum(subTabt.groupN == gi);
        
        caseINDS = find(subTabt.groupN == gi);
        
        %         % Thalamus Volume
        %         volDif = nan(numCASES,1);
        
        % Normalized Volume
        absDif = nan(numCASES,1);
        
        % Up or Down in size
        signDir = nan(numCASES,1);
        
        
        ti = 1;
        
        for si = 1:numCASES
            
            cI = caseINDS(si);
            fsurgC = subTabt.f_surg_n(cI);
            ssurgC = subTabt.s_surg_n(cI);
            
            
            
            fsurgInd = ismember(totalSegTable.Case,fsurgC);
            ssurgInd = ismember(totalSegTable.Case,ssurgC);
            
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
            
            
            
            if strcmp(sideIND,'L')
                
                FSindS = fsurgTab.Volmm3(ismember(fsurgTab.HemiS,'L') & ismember(fsurgTab.HippoArea,tmpBr))/FSICV;
                SSindS = ssurgTab.Volmm3(ismember(ssurgTab.HemiS,'L') & ismember(ssurgTab.HippoArea,tmpBr))/SSICV;
                
                FSindN = fsurgTab.Volmm3(ismember(fsurgTab.HemiS,'R') & ismember(fsurgTab.HippoArea,tmpBr))/FSICV;
                SSindN = ssurgTab.Volmm3(ismember(ssurgTab.HemiS,'R') & ismember(ssurgTab.HippoArea,tmpBr))/SSICV;
                
            else
                
                FSindS = fsurgTab.Volmm3(ismember(fsurgTab.HemiS,'R') & ismember(fsurgTab.HippoArea,tmpBr))/FSICV;
                SSindS = ssurgTab.Volmm3(ismember(ssurgTab.HemiS,'R') & ismember(ssurgTab.HippoArea,tmpBr))/SSICV;
                
                FSindN = fsurgTab.Volmm3(ismember(fsurgTab.HemiS,'L') & ismember(fsurgTab.HippoArea,tmpBr))/FSICV;
                SSindN = ssurgTab.Volmm3(ismember(ssurgTab.HemiS,'L') & ismember(ssurgTab.HippoArea,tmpBr))/SSICV;
                
            end
            
            if isempty(FSindS) || isempty(SSindS) || isempty(FSindN) || isempty(SSindN)
                groupIDS{gi,si} = nan;
                continue
            else
                
                groupIDS{gi,si} = fsurgC;
                FSdiff = (FSindS - FSindN) / max([FSindS , FSindN]);
                SSdiff = (SSindS - SSindN) / max([SSindS , SSindN]);
                
                absDif(ti,1) = FSdiff - SSdiff;
                signDir(ti,1) = (FSdiff - SSdiff)*100;
                
                
                ti = ti + 1;
            end
        end
        
        absDif = absDif(~isnan(absDif));
        signDir = signDir(~isnan(signDir));
        %         signDir = signDir(~isnan(signDir));
        
        dataOUTsub{1,gi} = absDif;
        dataOUTsub{2,gi} = signDir;
        %         dataOUTsub{3,gi} = signDir;
        
    end
    if outFlag == 1
        dataOUTmain{aai} = dataOUTsub(1,:);
    elseif outFlag == 2
        dataOUTmain{aai} = dataOUTsub(2,:);
    end
    
    
end
end











function [statsOUT] = getSTATS(allCaseDpd, groupNum, brNames)

statsOUT = cell(length(allCaseDpd),1);

for si = 1:length(allCaseDpd)
    
    forSTATS.data = [];
    forSTATS.group = [];
    forSTATS.brN = brNames{si};
    
    for g2 = 1:groupNum
        
        forSTATS.data = [forSTATS.data ; allCaseDpd{si,1}{1,g2}];
        
        grTmp = repmat(g2 , length(allCaseDpd{si,1}{1,g2}) , 1);
        
        forSTATS.group = [forSTATS.group ; grTmp];
        
    end
    
    statsOUT{si,1} = forSTATS;
    
end



end

