function [allCaseD , allSTATS] = subXFStabHemi_v3_HY()




bNames = {'entorhinal','lateralorbitofrontal','medialorbitofrontal'};

cd('Z:\Yilma_Project\CompiledCSVdata')

subTab = readtable('ob_subj_data.csv');

[allCaseD.PD] = getDATA(2, subTab, bNames, 'PD');
[allSTATS.PD] = getSTATS(allCaseD.PD, 2, bNames);



end




function [dataOUTmain] = getDATA(groupNum, subTab, brFlag, condI)

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
            
            if isempty(fsurgL) || isempty(fsurgR) || isempty(ssurgL) || isempty(ssurgR)
                continue
            else
                
                if strcmp(sideIND,'L')
                    
                    FSindS = fsurgL.SurfArea{ismember(fsurgL.StructName,tmpBr)};
                    SSindS = ssurgL.SurfArea{ismember(ssurgL.StructName,tmpBr)};
                    
                    FSindN = fsurgR.SurfArea{ismember(fsurgR.StructName,tmpBr)};
                    SSindN = ssurgR.SurfArea{ismember(ssurgR.StructName,tmpBr)};
                    
                else
                    
                    FSindS = fsurgR.SurfArea{ismember(fsurgR.StructName,tmpBr)};
                    SSindS = ssurgR.SurfArea{ismember(ssurgR.StructName,tmpBr)};
                    
                    FSindN = fsurgL.SurfArea{ismember(fsurgL.StructName,tmpBr)};
                    SSindN = ssurgL.SurfArea{ismember(ssurgL.StructName,tmpBr)};
                    
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

