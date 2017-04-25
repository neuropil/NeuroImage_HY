%% Step 1

% fprintf('Running Step 1 \n')
%
% FreeSurf_stats2txt_HY()
%
% fprintf('Step 1 Complete \n')

%% Step 2

% fprintf('Running Step 2 \n')
%
% FreeSurf_WMCombine_HY()
%
% fprintf('Step 2 Complete \n')

%% Step 3

% fprintf('Running Step 3 \n')
%
% CAT12_Extract_HY()
%
% fprintf('Step 3 Complete \n')

%% Step 4

function [outSTAT1pd, outSTAT2pd, outSTAT1etpd, outSTAT2etpd] = runHYgraphsStats_ab(ID)

switch ID
    case 'GM'
        
        [~ , inSTATS] = subXcatTab_v2_HY('GMmni');
        
    case 'WM'
        
        [~ , inSTATS] = subXcatTab_v2_HY('WMmni');
        
    case 'IH'
        
        [~ , inSTATS] = subXFStabCortM_v2_HY();
        
    case 'Caudate'
        
        [~ , inSTATS] = subXFStabSEG_v2_HY('Caud');
        
    case 'Hippocampus'
        
        [~ , inSTATS] = subXFStabSEG_v2_HY('Hippo');
        
    case 'Thalamus'
        
        [~ , inSTATS] = subXFStabSEG_v2_HY('Thal');

    case 'OBFc'
        
        [~ , inSTATS] = subXFStabHemi_v2_HY('latOB');

    case 'OB'
        
        [~ , inSTATS] = subOBVOL_v2_HY(1);

end


[~,~,outSTAT1pd] = kruskalwallis(inSTATS.PD{2}.data , inSTATS.PD{2}.group);
[outSTAT2pd,~,~] = multcompare(outSTAT1pd);
close all
boxplot(inSTATS.PD{2}.data , inSTATS.PD{2}.group)

set(gca,'YLim',[-0.1 4])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])

dataA = [inSTATS.PD{2}.data(inSTATS.PD{2}.group == 1) ; inSTATS.ET{2}.data(inSTATS.ET{2}.group == 1)];

groupA =[repmat({'p'},length(inSTATS.PD{2}.data(inSTATS.PD{2}.group == 1)),1) ;...
    repmat({'c'},length(inSTATS.ET{2}.group(inSTATS.ET{2}.group == 1)),1)];

[~,~,outSTAT1etpd] = ttest2(inSTATS.PD{2}.data(inSTATS.PD{2}.group == 1) , inSTATS.ET{2}.data(inSTATS.ET{2}.group == 1));
[outSTAT2etpd,~,~] = multcompare(outSTAT1etpd);
close all
boxplot(dataA , groupA)

set(gca,'YLim',[-0.05 5])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])


%% Group 3 long stimulation sub group analysis

%
% [GMallCaseD , GMallSTATS] = subXcatTab_v2_HY('GMmni');
% [WMallCaseD , WMallSTATS] = subXcatTab_v2_HY('WMmni');
% [impHemiDATA , impHemiSTATS] = subXFStabCortM_v2_HY();
%
% wm = WMallCaseD.PD{1,3};
% gm = GMallCaseD.PD{1,3};
% imp = impHemiDATA.PD{2,3};
%
% scatter(ones(size(gm)),gm); hold on
% scatter(ones(size(wm))+1,wm)
% scatter(ones(size(imp))+2,imp)
%
% xlim([0 4])

