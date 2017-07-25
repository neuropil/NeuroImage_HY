function [allCaseD , allSTATS] = subXhippo_HY()




cd('Z:\Yilma_Project\CompiledCSVdata')

subTab = readtable('ob_subj_data.csv');

load('TotalHippoTable.mat');

hippoSegTable.Case = cellfun(@(x) str2double(x), hippoSegTable.Case); %#ok<NODEF>

bNames = unique(hippoSegTable.HippoArea);

[allCaseD.PD] = getDATA(2, subTab, hippoSegTable, bNames, 'PD');
[allSTATS.PD] = getSTATS(allCaseD.PD, 2 , bNames);




end






function [dataOUTmain] = getDATA(groupNum, subTab, totalSegTable, brFlag, condI)

condIn = ismember(subTab.cond,condI);
subTabt = subTab(condIn,:);

dataOUTmain = cell(length(brFlag),1);
for aai = 1:length(brFlag)
    
    tmpBr = brFlag{aai};
    
    dataOUTsub = cell(2,groupNum);
    
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
            
            if strcmp(sideIND,'L')
                
                FSindS = fsurgTab.Volmm3(ismember(fsurgTab.HemiS,'L') & ismember(fsurgTab.HippoArea,tmpBr));
                SSindS = ssurgTab.Volmm3(ismember(ssurgTab.HemiS,'L') & ismember(ssurgTab.HippoArea,tmpBr));
                
                FSindN = fsurgTab.Volmm3(ismember(fsurgTab.HemiS,'R') & ismember(fsurgTab.HippoArea,tmpBr));
                SSindN = ssurgTab.Volmm3(ismember(ssurgTab.HemiS,'R') & ismember(ssurgTab.HippoArea,tmpBr));
                
            else
                
                FSindS = fsurgTab.Volmm3(ismember(fsurgTab.HemiS,'R') & ismember(fsurgTab.HippoArea,tmpBr));
                SSindS = ssurgTab.Volmm3(ismember(ssurgTab.HemiS,'R') & ismember(ssurgTab.HippoArea,tmpBr));
                
                FSindN = fsurgTab.Volmm3(ismember(fsurgTab.HemiS,'L') & ismember(fsurgTab.HippoArea,tmpBr));
                SSindN = ssurgTab.Volmm3(ismember(ssurgTab.HemiS,'L') & ismember(ssurgTab.HippoArea,tmpBr));
                
            end
            
            if isempty(FSindS) || isempty(SSindS) || isempty(FSindN) || isempty(SSindN)
                continue
            else
                
%                 volDif(ti,1) = double(ssurgTab.Volmm3(SSindT) - fsurgTab.Volmm3(FSindT)) / double(fsurgTab.Volmm3(FSindT));
                
                FSdiff = (FSindS - FSindN) / max([FSindS , FSindN]);
                SSdiff = (SSindS - SSindN) / max([SSindS , SSindN]);
                
                absDif(ti,1) = abs(FSdiff - SSdiff)*100;
                signDir(ti,1) = sign(FSdiff - SSdiff);
                 
                
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
        
        forSTATS.data = [forSTATS.data ; allCaseDpd{si,1}{1,g2} .* allCaseDpd{si,1}{2,g2}];
        
        grTmp = repmat(g2 , length(allCaseDpd{si,1}{1,g2}) , 1);
        
        forSTATS.group = [forSTATS.group ; grTmp];
        
    end
    
    statsOUT{si,1} = forSTATS;
    
end



end

