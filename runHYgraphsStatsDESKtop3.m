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


function [outSTAT1pd, outSTAT2pd, outSTAT1etpd] = runHYgraphsStatsDESKtop3(ID)

switch ID
    case 'GM' %%% DONE
        
        [~ , GMallSTATS] = subXcatTab_v3_HY('GMmni');
        
        [outSTAT1pd,~,~] = ranksum(GMallSTATS.PD{1}.data(GMallSTATS.PD{1}.group == 1) , GMallSTATS.PD{1}.data(GMallSTATS.PD{1}.group == 2));
        
        close all
        boxplot(GMallSTATS.PD{1}.data , GMallSTATS.PD{1}.group)
        
%         set(gca,'YLim',[0 0.05])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
        
    case 'WM' %%% DONE
        
        [~ , WMallSTATS] = subXcatTab_v3_HY('WMmni');
        
        [outSTAT1pd,~,~] = ranksum(WMallSTATS.PD{1}.data(WMallSTATS.PD{1}.group == 1) , WMallSTATS.PD{1}.data(WMallSTATS.PD{1}.group == 2));
        
        close all
        boxplot(WMallSTATS.PD{1}.data , WMallSTATS.PD{1}.group)
        
%         set(gca,'YLim',[0 0.09])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
    case 'IH' %%% DONE
        
        [~ , impHemiSTATS] = subXFStabCortM_v3_HY();
        
        [outSTAT1pd,~,~] = ranksum(impHemiSTATS.PD{1}.data(impHemiSTATS.PD{1}.group == 1) , impHemiSTATS.PD{1}.data(impHemiSTATS.PD{1}.group == 2));
        
        close all
        boxplot(impHemiSTATS.PD{1}.data , impHemiSTATS.PD{1}.group)
        
        set(gca,'YLim',[-0.001 0.16])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
        
        
    case 'Hippocampus' %%% DONE
        
        [~ , hidallSTATS] = subXhippo_HY();
        
        for pi = 1:length(hidallSTATS.PD)
            figure
            
            tmpData = hidallSTATS.PD{pi}.data;
            tmpID = hidallSTATS.PD{pi}.group;
            
            [pval,~,~ ] = ranksum(tmpData(tmpID == 1),tmpData(tmpID == 2));
            
            boxplot(tmpData,tmpID)
            
            set(gca,'XTick',1:2)
            set(gca,'XTickLabel',{'NoStim','Stim'})
            title([hidallSTATS.PD{pi}.brN ,'  ', num2str(pval)])
            
        end
        
    case 'Thal-Caud' %%% DONE
        
        [~ , thdallSTATS] = subXFStabSEG_v3_HY();
        close all
        for pi = 1:length(thdallSTATS.PD)
            figure
            
            tmpData = thdallSTATS.PD{pi}.data;
            tmpID = thdallSTATS.PD{pi}.group;
            
            [pval,~,~ ] = ranksum(tmpData(tmpID == 1),tmpData(tmpID == 2));
            
            boxplot(tmpData,tmpID)
            
            set(gca,'XTick',1:2)
            set(gca,'XTickLabel',{'NoStim','Stim'})
            title([thdallSTATS.PD{pi}.brN ,'  ', num2str(pval)])
            
        end
        
        
    case 'OBFc' %%% DONE
        
        [~ , loballSTATS] = subXFStabHemi_v3_HY();
        
        close all
        for pi = 1:length(loballSTATS.PD)
            figure
            
            tmpData = loballSTATS.PD{pi}.data;
            tmpID = loballSTATS.PD{pi}.group;
            
            [pval,~,~ ] = ranksum(tmpData(tmpID == 1),tmpData(tmpID == 2));
            
            boxplot(tmpData,tmpID)
            
            set(gca,'XTick',1:2)
            set(gca,'XTickLabel',{'NoStim','Stim'})
            title([loballSTATS.PD{pi}.brN ,'  ', num2str(pval)])
            
        end
        
        
    case 'OB' %%% DONE
        
        %         [~ , OBallSTATS] = subOBVOL_v2_HY(1);
        
        [~ , OBallSTATS] = subOBVOL_v2b_HY(1);
        
        [outSTAT2pd,~,~] = ranksum(OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 1) , OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 2));
        
        close all
        boxplot(OBallSTATS.PD{2}.data , OBallSTATS.PD{2}.group)
        
        set(gca,'XTick',1:2)
        set(gca,'XTickLabel',{'NoStim','Stim'})
        
        
        set(gca,'YLim',[-0.1 0.1])
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
        g1 = OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 1);
        g2 = OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 2);
        ga1 = abs(g1);
        ga2 = abs(g2);
        t = readtable('ob_subj_data.csv');
        tPD = t(ismember(t.cond,'PD'),:);
        tPDg2 = tPD(tPD.groupN == 2,:);
        tPDg1 = tPD(tPD.groupN == 1,:);
        stD = tPDg2.stimdur;
        nstD = tPDg1.stimdur;
        
        figure;
        plot(stD,g2,'ro')
        hold on
        line([min(get(gca,'XTick')) , max(get(gca,'XTick'))],[0 0],'LineStyle','--','Color','k')
        lsline
        [rho , pval] = corr(stD,g2);
        
        figure;
        plot(stD,ga2,'ro')
        hold on
        line([min(get(gca,'XTick')) , max(get(gca,'XTick'))],[0 0],'LineStyle','--','Color','k')
        lsline
        plot(nstD,ga1,'bo')
        [rho , pval] = corr(stD,ga2);
        
        figure;
        plot(stD,g2,'ro')
        hold on
        line([min(get(gca,'XTick')) , max(get(gca,'XTick'))],[0 0],'LineStyle','--','Color','k')
        lsline
        [rho , pval] = corr(stD,g2);
        
        
        
        
        
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

