function [allCaseD , allSTATS] = subOBVOL_v2b_HY(absFlag)

cd('Z:\Yilma_Project\CompiledCSVdata')

load('AllOBdata.mat');

allOBTab.caseID = cellfun(@(x) str2double(x(2:4)), allOBTab.caseID); %#ok<NODEF>

cd('Z:\Yilma_Project\CompiledCSVdata')

subTab = readtable('ob_subj_data.csv');



            [allCaseD.PD] = getDATA(2, subTab, allOBTab, absFlag, 'PD');
            [allSTATS.PD] = getSTATS(allCaseD.PD, 2);





end


function [dataOUT] = getDATA(groupNum, subTab, allOBTab, absFlag, condI)

condIn = ismember(subTab.cond,condI);
subTabt = subTab(condIn,:);

dataOUT = cell(1,groupNum);

for gi = 1:groupNum
    
    numCASES = sum(subTabt.groupN == gi);
    
    caseINDS = find(subTabt.groupN == gi);
    
    % Thalamus Volume
    volDif = nan(numCASES,1);
    
    % Normalized Volume
    normDif = nan(numCASES,1);
    
    ti = 1;
    
    for si = 1:numCASES
        
        cI = caseINDS(si);
        fsurgC = subTabt.f_surg_n(cI);
        ssurgC = subTabt.s_surg_n(cI);
        
        fsurgInd = ismember(allOBTab.caseID,fsurgC);
        ssurgInd = ismember(allOBTab.caseID,ssurgC);
        
        fsurgTab = allOBTab(fsurgInd,:);
        ssurgTab = allOBTab(ssurgInd,:);
        
        sideIND = subTabt.f_surg_s{cI};
        
        if strcmp(sideIND,'L')
            stimFlag = 'L_OB';
            nstimFlag = 'R_OB';
        else
            stimFlag = 'R_OB';
            nstimFlag = 'L_OB';
        end
        
        if sum(fsurgInd) == 0 || sum(ssurgInd) == 0
            continue
        else
            
            FSIndS = fsurgTab.Volume_mm3(ismember(fsurgTab.LabelName,stimFlag));
            SSIndS = ssurgTab.Volume_mm3(ismember(ssurgTab.LabelName,stimFlag));
            
            FSindN = fsurgTab.Volume_mm3(ismember(fsurgTab.LabelName,nstimFlag));
            SSIndN = ssurgTab.Volume_mm3(ismember(ssurgTab.LabelName,nstimFlag));
            
            if absFlag == 1
                
                totalFS = sum(double(fsurgTab.Volume_mm3(:)));
                fracFS = FSIndS/totalFS;
                
                totalSS = sum(double(ssurgTab.Volume_mm3(:)));
                fracSS = SSIndS/totalSS;
                
                diffPrePost = fracSS - fracFS;
                
                volDif(ti,1) = diffPrePost;
                
                
                FsDiff = (FSIndS - FSindN) / max([FSindN , FSIndS]);
                SsDiff = (SSIndS - SSIndN) / max([SSIndS , SSIndN]);
                
                perDif = (SsDiff - FsDiff)*100;
                
                
            elseif absFlag == 2
                volDif(ti,1) = abs(double(ssurgTab.Volume_mm3(SSindT) - fsurgTab.Volume_mm3(FSindT))) / double(fsurgTab.Volume_mm3(FSindT));
            elseif absFlag == 3
                volDif(ti,1) = double(ssurgTab.Volume_mm3(SSindT) - fsurgTab.Volume_mm3(FSindT)) / double(fsurgTab.Volume_mm3(FSindT));
            end
            %             normDif(ti,1) = abs(double(ssurgTab.IMmean(SSindT) - fsurgTab.IMmean(FSindT))) / double(fsurgTab.IMmean(FSindT));
            normDif(ti,1) = perDif;
            ti = ti + 1;
        end
    end
    
    volDif = volDif(~isnan(volDif));
    normDif = normDif(~isnan(normDif));
    
    
    dataOUT{1,gi} = volDif;
    dataOUT{2,gi} = normDif;
    
end


end



function [statsOUT] = getSTATS(allCaseDpd, groupNum)

statsOUT = cell(1,groupNum);

for si = 1:2
    
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


