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

% OB_VolumeExtract_HY()


%% Step 5

% OB_CreateAllCSV_HY()

%% Step 6

% This function


function [outSTAT1pd, outSTAT2pd, outSTAT1etpd] = runHYgraphsStats(ID)

switch ID
    case 'GM'
        
        [~ , GMallSTATS] = subXcatTab_v2_HY('GMmni');
        
        [~,~,outSTAT1pd] = kruskalwallis(GMallSTATS.PD{1}.data , GMallSTATS.PD{1}.group);
        [outSTAT2pd,~,~] = multcompare(outSTAT1pd);
        close all
        boxplot(GMallSTATS.PD{1}.data , GMallSTATS.PD{1}.group)
        
        set(gca,'YLim',[0 0.05])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
        dataA = [GMallSTATS.PD{1}.data(GMallSTATS.PD{1}.group == 1) ; GMallSTATS.ET{1}.data(GMallSTATS.ET{1}.group == 1)];
        
        groupA =[repmat({'p'},length(GMallSTATS.PD{1}.data(GMallSTATS.PD{1}.group == 1)),1) ;...
            repmat({'c'},length(GMallSTATS.ET{1}.group(GMallSTATS.ET{1}.group == 1)),1)];
        
        [outSTAT1etpd,~,~] = ranksum(GMallSTATS.PD{1}.data(GMallSTATS.PD{1}.group == 1) , GMallSTATS.ET{1}.data(GMallSTATS.ET{1}.group == 1));
        figure;
        boxplot(dataA , groupA)
        
        set(gca,'YLim',[0 0.05])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
        
    case 'WM'
        
        [~ , WMallSTATS] = subXcatTab_v2_HY('WMmni');
        
        [~,~,outSTAT1pd] = kruskalwallis(WMallSTATS.PD{1}.data , WMallSTATS.PD{1}.group);
        [outSTAT2pd,~,~] = multcompare(outSTAT1pd);
        close all
        boxplot(WMallSTATS.PD{1}.data , WMallSTATS.PD{1}.group)
        
        set(gca,'YLim',[0 0.09])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
        dataA = [WMallSTATS.PD{1}.data(WMallSTATS.PD{1}.group == 1) ; WMallSTATS.ET{1}.data(WMallSTATS.ET{1}.group == 1)];
        
        groupA =[repmat({'p'},length(WMallSTATS.PD{1}.data(WMallSTATS.PD{1}.group == 1)),1) ;...
            repmat({'c'},length(WMallSTATS.ET{1}.group(WMallSTATS.ET{1}.group == 1)),1)];
        
        [outSTAT1etpd,~,~] = ranksum(WMallSTATS.PD{1}.data(WMallSTATS.PD{1}.group == 1) , WMallSTATS.ET{1}.data(WMallSTATS.ET{1}.group == 1));
        %         [outSTAT2etpd,~,~] = multcompare(outSTAT1etpd);
        figure;
        boxplot(dataA , groupA)
        
        set(gca,'YLim',[0 0.09])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
    case 'IH'
        
        [~ , impHemiSTATS] = subXFStabCortM_v2_HY();
        
        [~,~,outSTAT1pd] = kruskalwallis(impHemiSTATS.PD{2}.data , impHemiSTATS.PD{2}.group);
        [outSTAT2pd,~,~] = multcompare(outSTAT1pd);
        close all
        boxplot(impHemiSTATS.PD{2}.data , impHemiSTATS.PD{2}.group)
        
        set(gca,'YLim',[-0.001 0.16])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
        dataA = [impHemiSTATS.PD{2}.data(impHemiSTATS.PD{2}.group == 1) ; impHemiSTATS.ET{2}.data(impHemiSTATS.ET{2}.group == 1)];
        
        groupA =[repmat({'p'},length(impHemiSTATS.PD{2}.data(impHemiSTATS.PD{2}.group == 1)),1) ;...
            repmat({'c'},length(impHemiSTATS.ET{2}.group(impHemiSTATS.ET{2}.group == 1)),1)];
        
        [outSTAT1etpd,~,~] = ranksum(impHemiSTATS.PD{2}.data(impHemiSTATS.PD{2}.group == 1) , impHemiSTATS.ET{2}.data(impHemiSTATS.ET{2}.group == 1));
        %         [outSTAT2etpd,~,~] = multcompare(outSTAT1etpd);
        figure;
        boxplot(dataA , groupA)
        
        set(gca,'YLim',[-0.001 0.16])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
    case 'Caudate'
        [~ , cadallSTATS] = subXFStabSEG_v2_HY('Caud');
        
        [~,~,outSTAT1pd] = kruskalwallis(cadallSTATS.PD{2}.data , cadallSTATS.PD{2}.group);
        [outSTAT2pd,~,~] = multcompare(outSTAT1pd);
        close all
        boxplot(cadallSTATS.PD{2}.data , cadallSTATS.PD{2}.group)
        
        set(gca,'YLim',[-0.001 0.16])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
        dataA = [cadallSTATS.PD{2}.data(cadallSTATS.PD{2}.group == 1) ; cadallSTATS.ET{2}.data(cadallSTATS.ET{2}.group == 1)];
        
        groupA =[repmat({'p'},length(cadallSTATS.PD{2}.data(cadallSTATS.PD{2}.group == 1)),1) ;...
            repmat({'c'},length(cadallSTATS.ET{2}.group(cadallSTATS.ET{2}.group == 1)),1)];
        
        [outSTAT1etpd,~,~] = ranksum(cadallSTATS.PD{2}.data(cadallSTATS.PD{2}.group == 1) , cadallSTATS.ET{2}.data(cadallSTATS.ET{2}.group == 1));
        
        figure;
        boxplot(dataA , groupA)
        
        set(gca,'YLim',[-0.001 0.16])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
    case 'Hippocampus'
        
        [~ , hidallSTATS] = subXFStabSEG_v2_HY('Hippo');
        
        [~,~,outSTAT1pd] = kruskalwallis(hidallSTATS.PD{2}.data , hidallSTATS.PD{2}.group);
        [outSTAT2pd,~,~] = multcompare(outSTAT1pd);
        close all
        boxplot(hidallSTATS.PD{2}.data , hidallSTATS.PD{2}.group)
        
        set(gca,'YLim',[-0.001 0.3])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
        dataA = [hidallSTATS.PD{2}.data(hidallSTATS.PD{2}.group == 1) ; hidallSTATS.ET{2}.data(hidallSTATS.ET{2}.group == 1)];
        
        groupA =[repmat({'p'},length(hidallSTATS.PD{2}.data(hidallSTATS.PD{2}.group == 1)),1) ;...
            repmat({'c'},length(hidallSTATS.ET{2}.group(hidallSTATS.ET{2}.group == 1)),1)];
        
        [outSTAT1etpd,~,~] = ranksum(hidallSTATS.PD{2}.data(hidallSTATS.PD{2}.group == 1) , hidallSTATS.ET{2}.data(hidallSTATS.ET{2}.group == 1));
        
        figure;
        boxplot(dataA , groupA)
        
        set(gca,'YLim',[-0.001 0.3])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
    case 'Thalamus'
        
        [~ , thdallSTATS] = subXFStabSEG_v2_HY('Thal');
        
        [~,~,outSTAT1pd] = kruskalwallis(thdallSTATS.PD{2}.data , thdallSTATS.PD{2}.group);
        [outSTAT2pd,~,~] = multcompare(outSTAT1pd);
        close all
        boxplot(thdallSTATS.PD{2}.data , thdallSTATS.PD{2}.group)
        
        set(gca,'YLim',[-0.001 0.2])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
        dataA = [thdallSTATS.PD{2}.data(thdallSTATS.PD{2}.group == 1) ; thdallSTATS.ET{2}.data(thdallSTATS.ET{2}.group == 1)];
        
        groupA =[repmat({'p'},length(thdallSTATS.PD{2}.data(thdallSTATS.PD{2}.group == 1)),1) ;...
            repmat({'c'},length(thdallSTATS.ET{2}.group(thdallSTATS.ET{2}.group == 1)),1)];
        
        [outSTAT1etpd,~,~] = ranksum(thdallSTATS.PD{2}.data(thdallSTATS.PD{2}.group == 1) , thdallSTATS.ET{2}.data(thdallSTATS.ET{2}.group == 1));

        figure;
        boxplot(dataA , groupA)
        
        set(gca,'YLim',[-0.001 0.2])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
    case 'OBFc'
        
        [~ , loballSTATS] = subXFStabHemi_v2_HY('latOB');
        
        [~,~,outSTAT1pd] = kruskalwallis(loballSTATS.PD{1}.data , loballSTATS.PD{1}.group);
        [outSTAT2pd,~,~] = multcompare(outSTAT1pd);
        close all
        boxplot(loballSTATS.PD{1}.data , loballSTATS.PD{1}.group)
        
        set(gca,'YLim',[-0.001 0.3])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
        dataA = [loballSTATS.PD{1}.data(loballSTATS.PD{1}.group == 1) ; loballSTATS.ET{1}.data(loballSTATS.ET{1}.group == 1)];
        
        groupA =[repmat({'p'},length(loballSTATS.PD{1}.data(loballSTATS.PD{1}.group == 1)),1) ;...
            repmat({'c'},length(loballSTATS.ET{1}.group(loballSTATS.ET{1}.group == 1)),1)];
        
        [~,outSTAT1etpd,~] = ttest2(loballSTATS.PD{1}.data(loballSTATS.PD{1}.group == 1) , loballSTATS.ET{1}.data(loballSTATS.ET{1}.group == 1));

        figure;
        boxplot(dataA , groupA)
        
        set(gca,'YLim',[-0.001 0.3])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
    case 'OB'
        
%         [~ , OBallSTATS] = subOBVOL_v2_HY(1);
        
        [~ , OBallSTATS] = subOBVOL_v2b_HY(1);
        [~,~,outSTAT1pd] = kruskalwallis(OBallSTATS.PD{1}.data , OBallSTATS.PD{1}.group);
        [outSTAT2pd,~,~] = multcompare(outSTAT1pd);
        close all
        boxplot(OBallSTATS.PD{1}.data , OBallSTATS.PD{1}.group)
        
        set(gca,'YLim',[-0.1 0.1])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
        dataA = [OBallSTATS.PD{1}.data(OBallSTATS.PD{1}.group == 1) ; OBallSTATS.ET{1}.data(OBallSTATS.ET{1}.group == 1)];
        
        groupA =[repmat({'p'},length(OBallSTATS.PD{1}.data(OBallSTATS.PD{1}.group == 1)),1) ;...
            repmat({'c'},length(OBallSTATS.ET{2}.group(OBallSTATS.ET{1}.group == 1)),1)];
        
        [~,outSTAT1etpd,~] = ttest2(OBallSTATS.PD{1}.data(OBallSTATS.PD{1}.group == 1) , OBallSTATS.ET{1}.data(OBallSTATS.ET{1}.group == 1));

        figure;
        boxplot(dataA , groupA)
        
        set(gca,'YLim',[-0.1 0.1])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
        
        
end





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

