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

function

[GMallCaseD , GMallSTATS] = subXcatTab_v2_HY('GMmni');

[a,b,c] = kruskalwallis(GMallSTATS.PD{1}.data , GMallSTATS.PD{1}.group);
[a,b,c] = multcompare(c)
close all
boxplot(GMallSTATS.PD{1}.data , GMallSTATS.PD{1}.group)

set(gca,'YLim',[-0.005 0.05])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])



%%

[WMallCaseD , WMallSTATS] = subXcatTab_v2_HY('WMmni');

[a,b,c] = kruskalwallis(WMallSTATS.PD{1}.data , WMallSTATS.PD{1}.group);
[a,b,c] = multcompare(c)
close all
boxplot(WMallSTATS.PD{1}.data , WMallSTATS.PD{1}.group)

set(gca,'YLim',[-0.005 0.09])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])


%% GM - PD vs ET group 1 

dataA = [GMallSTATS.PD{1}.data(GMallSTATS.PD{1}.group == 1) ; GMallSTATS.ET{1}.data(GMallSTATS.ET{1}.group == 1)]

groupA =[repmat({'p'},length(GMallSTATS.PD{1}.data(GMallSTATS.PD{1}.group == 1)),1) ;...
    repmat({'c'},length(GMallSTATS.ET{1}.group(GMallSTATS.ET{1}.group == 1)),1)]

[a,b,c] = ranksum(GMallSTATS.PD{1}.data(GMallSTATS.PD{1}.group == 1) , GMallSTATS.ET{1}.data(GMallSTATS.ET{1}.group == 1));
[a,b,c] = multcompare(c)
close all
boxplot(dataA , groupA)

set(gca,'YLim',[-0.005 0.09])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])

%% WM - PD vs ET group 1 

dataA = [WMallSTATS.PD{1}.data(WMallSTATS.PD{1}.group == 1) ; WMallSTATS.ET{1}.data(WMallSTATS.ET{1}.group == 1)]

groupA =[repmat({'p'},length(WMallSTATS.PD{1}.data(WMallSTATS.PD{1}.group == 1)),1) ;...
    repmat({'c'},length(WMallSTATS.ET{1}.group(WMallSTATS.ET{1}.group == 1)),1)]

[a,b,c] = ranksum(WMallSTATS.PD{1}.data(WMallSTATS.PD{1}.group == 1) , WMallSTATS.ET{1}.data(WMallSTATS.ET{1}.group == 1));

close all
boxplot(dataA , groupA)

set(gca,'YLim',[-0.005 0.09])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])


%% Implanted Hemisphere

[impHemiDATA , impHemiSTATS] = subXFStabCortM_v2_HY()


[a,b,c] = kruskalwallis(impHemiSTATS.PD{2}.data , impHemiSTATS.PD{2}.group);
[a,b,c] = multcompare(c)
close all
boxplot(impHemiSTATS.PD{2}.data , impHemiSTATS.PD{2}.group)

set(gca,'YLim',[-0.005 0.16])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])

%%

dataA = [impHemiSTATS.PD{2}.data(impHemiSTATS.PD{2}.group == 1) ; impHemiSTATS.ET{2}.data(impHemiSTATS.ET{2}.group == 1)]

groupA =[repmat({'p'},length(impHemiSTATS.PD{2}.data(impHemiSTATS.PD{2}.group == 1)),1) ;...
    repmat({'c'},length(impHemiSTATS.ET{2}.group(impHemiSTATS.ET{2}.group == 1)),1)]

[a,b,c] = ranksum(impHemiSTATS.PD{2}.data(impHemiSTATS.PD{2}.group == 1) , impHemiSTATS.ET{2}.data(impHemiSTATS.ET{2}.group == 1))

close all
boxplot(dataA , groupA)

set(gca,'YLim',[-0.005 0.09])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])

%% PD Caudate Volume 

[cadallCaseD , cadallSTATS] = subXFStabSEG_v2_HY('Caud')

[a,b,c] = kruskalwallis(cadallSTATS.PD{2}.data , cadallSTATS.PD{2}.group);
[a,b,c] = multcompare(c)
close all
boxplot(cadallSTATS.PD{2}.data , cadallSTATS.PD{2}.group)

set(gca,'YLim',[-0.005 0.16])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])

%% PD Caudate Volume vs ET

dataA = [cadallSTATS.PD{2}.data(cadallSTATS.PD{2}.group == 1) ; cadallSTATS.ET{2}.data(cadallSTATS.ET{2}.group == 1)]

groupA =[repmat({'p'},length(cadallSTATS.PD{2}.data(cadallSTATS.PD{2}.group == 1)),1) ;...
    repmat({'c'},length(cadallSTATS.ET{2}.group(cadallSTATS.ET{2}.group == 1)),1)]

[a,b,c] = ranksum(cadallSTATS.PD{2}.data(cadallSTATS.PD{2}.group == 1) , cadallSTATS.ET{2}.data(cadallSTATS.ET{2}.group == 1))

close all
boxplot(dataA , groupA)

set(gca,'YLim',[-0.005 0.16])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])

%% PD Hippocampus Volume 

[hidallCaseD , hidallSTATS] = subXFStabSEG_v2_HY('Hippo')

[a,b,c] = kruskalwallis(hidallSTATS.PD{2}.data , hidallSTATS.PD{2}.group);
[a,b,c] = multcompare(c)
close all
boxplot(hidallSTATS.PD{2}.data , hidallSTATS.PD{2}.group)

set(gca,'YLim',[-0.005 0.28])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])


%% PD Hippocampus Volume vs ET

dataA = [hidallSTATS.PD{2}.data(hidallSTATS.PD{2}.group == 1) ; hidallSTATS.ET{2}.data(hidallSTATS.ET{2}.group == 1)]

groupA =[repmat({'p'},length(hidallSTATS.PD{2}.data(hidallSTATS.PD{2}.group == 1)),1) ;...
    repmat({'c'},length(hidallSTATS.ET{2}.group(hidallSTATS.ET{2}.group == 1)),1)]

[a,b,c] = ranksum(hidallSTATS.PD{2}.data(hidallSTATS.PD{2}.group == 1) , hidallSTATS.ET{2}.data(hidallSTATS.ET{2}.group == 1))

close all
boxplot(dataA , groupA)

set(gca,'YLim',[-0.005 0.16])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])

%% PD thalamus Volume 

[thdallCaseD , thdallSTATS] = subXFStabSEG_v2_HY('Thal')

[a,b,c] = kruskalwallis(thdallSTATS.PD{2}.data , thdallSTATS.PD{2}.group);
[a,b,c] = multcompare(c)
close all
boxplot(thdallSTATS.PD{2}.data , thdallSTATS.PD{2}.group)

set(gca,'YLim',[-0.005 0.18])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])


%% PD thalamus Volume vs ET

dataA = [thdallSTATS.PD{2}.data(thdallSTATS.PD{2}.group == 1) ; thdallSTATS.ET{2}.data(thdallSTATS.ET{2}.group == 1)]

groupA =[repmat({'p'},length(thdallSTATS.PD{2}.data(thdallSTATS.PD{2}.group == 1)),1) ;...
    repmat({'c'},length(thdallSTATS.ET{2}.group(thdallSTATS.ET{2}.group == 1)),1)]

[a,b,c] = ranksum(thdallSTATS.PD{2}.data(thdallSTATS.PD{2}.group == 1) , thdallSTATS.ET{2}.data(thdallSTATS.ET{2}.group == 1))

close all
boxplot(dataA , groupA)

set(gca,'YLim',[-0.005 0.1])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])

%% PD lateral orbitofrontal 

[loballCaseD , loballSTATS] = subXFStabHemi_v2_HY('latOB')

[a,b,c] = kruskalwallis(loballSTATS.PD{1}.data , loballSTATS.PD{1}.group);
[a,b,c] = multcompare(c)
close all
boxplot(loballSTATS.PD{1}.data , loballSTATS.PD{1}.group)

set(gca,'YLim',[-0.005 0.28])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])


%% PD lateral orbitofrontal v ET

dataA = [loballSTATS.PD{1}.data(loballSTATS.PD{1}.group == 1) ; loballSTATS.ET{1}.data(loballSTATS.ET{1}.group == 1)]

groupA =[repmat({'p'},length(loballSTATS.PD{1}.data(loballSTATS.PD{1}.group == 1)),1) ;...
    repmat({'c'},length(loballSTATS.ET{1}.group(loballSTATS.ET{1}.group == 1)),1)]

[a,b,c] = ttest2(loballSTATS.PD{1}.data(loballSTATS.PD{1}.group == 1) , loballSTATS.ET{1}.data(loballSTATS.ET{1}.group == 1))

close all
boxplot(dataA , groupA)

set(gca,'YLim',[-0.005 0.15])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])

%%

[OBallCaseD , OBallSTATS] = subOBVOL_v2_HY(1);

[a,b,c] = kruskalwallis(OBallSTATS.PD{2}.data , OBallSTATS.PD{2}.group);
[a,b,c] = multcompare(c)
close all
boxplot(OBallSTATS.PD{2}.data , OBallSTATS.PD{2}.group)

set(gca,'YLim',[-0.1 4])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])


%% PD lateral orbitofrontal v ET

dataA = [OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 1) ; OBallSTATS.ET{2}.data(OBallSTATS.ET{2}.group == 1)]

groupA =[repmat({'p'},length(OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 1)),1) ;...
    repmat({'c'},length(OBallSTATS.ET{2}.group(OBallSTATS.ET{2}.group == 1)),1)]

[a,b,c] = ttest2(OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 1) , OBallSTATS.ET{2}.data(OBallSTATS.ET{2}.group == 1))

close all
boxplot(dataA , groupA)

set(gca,'YLim',[-0.05 5])
yVALs = get(gca,'YTick');

yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])









%% Group 3 long stimulation sub group analysis


[GMallCaseD , GMallSTATS] = subXcatTab_v2_HY('GMmni');
[WMallCaseD , WMallSTATS] = subXcatTab_v2_HY('WMmni');
[impHemiDATA , impHemiSTATS] = subXFStabCortM_v2_HY();

wm = WMallCaseD.PD{1,3};
gm = GMallCaseD.PD{1,3};
imp = impHemiDATA.PD{2,3};

scatter(ones(size(gm)),gm); hold on
scatter(ones(size(wm))+1,wm)
scatter(ones(size(imp))+2,imp)

xlim([0 4])

