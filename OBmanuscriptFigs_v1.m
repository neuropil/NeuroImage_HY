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


function [outSTAT1pd, outSTAT2pd, outSTAT1etpd] = OBmanuscriptFigs_v1(ID)

switch ID
    
    case 'Hippocampus' %%% DONE
        
        [~ , hidallSTATS, ids] = subXhippo_HYb(2,1);
        gr2IDS = transpose([ids{2,:}]);
        redS = linspace(0,0.9,length(hidallSTATS.PD));
        
        for pi = 1:length(hidallSTATS.PD)
            figure
            
            tmpData = hidallSTATS.PD{pi}.data;
            tmpID = hidallSTATS.PD{pi}.group;
            
            [pval,~,~ ] = ranksum(tmpData(tmpID == 1),tmpData(tmpID == 2));
            
            boxplot(tmpData,tmpID,'Colors',[0 0 0])
            
            yVALs = get(gca,'YTick');
            
            yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
            
            set(gca,'XTick',1:2)
            set(gca,'XTickLabel',{'NoStim','Stim'})
            title([hidallSTATS.PD{pi}.brN],'Interpreter','none')
            
            scatDat = [[ hidallSTATS.PD{1}.data(hidallSTATS.PD{1}.group == 1) ;...
                hidallSTATS.PD{1}.data(hidallSTATS.PD{1}.group == 2) ],...
                [ones(size(hidallSTATS.PD{1}.data(hidallSTATS.PD{1}.group == 1))) ;...
                repmat(2,size(hidallSTATS.PD{1}.data(hidallSTATS.PD{1}.group == 2)))]];
            
            scatCols = [repmat([0 0 0],size(hidallSTATS.PD{1}.data(hidallSTATS.PD{1}.group == 1))) ;...
                repmat([1 redS(pi) redS(pi)],size(hidallSTATS.PD{1}.data(hidallSTATS.PD{1}.group == 2)))];
            
            hold on
            
            scatter(scatDat(:,2),scatDat(:,1),40,scatCols,'filled')
            
            %% Add p-value
            text(2,max(yVALs) - abs(min(yVALs)/2), sprintf('p = %0.3f',pval))
            
            axis square
            
            cd('E:\Dropbox\Publications_Meta\InProgress\Yilma_OB_DBS\images\NewFigs08032017')
            
            fname = ['Figure1_',hidallSTATS.PD{pi}.brN,'.pdf'];
            
            saveas(gcf, fname)
            
        end
        close all
        figure
        
        tmpMin = 0;
        tmpMax = 5;
        for pi2 = 1:length(hidallSTATS.PD)
            
            cd('Z:\Yilma_Project\CompiledCSVdata')
            
            %         g1 = hidallSTATS.PD{2}.data(hidallSTATS.PD{2}.group == 1);
            g2 = hidallSTATS.PD{pi2}.data(hidallSTATS.PD{pi2}.group == 2);
            %         ga1 = abs(g1);
            ga2 = abs(g2);
            t = readtable('ob_subj_data.csv');
            tPD = t(ismember(t.cond,'PD'),:);
            tPDg2 = tPD(tPD.groupN == 2,:);
            %         tPDg1 = tPD(tPD.groupN == 1,:);
            
            tPDg3 = tPDg2(tPDg2.f_surg_n == gr2IDS,:);
            stD = tPDg3.stimdur;
            %         nstD = tPDg1.stimdur;
            title('hippocampus')
            coef_fit = polyfit(stD,ga2,1);
            y_fit = polyval(coef_fit,stD);
            
            line([min(stD) max(stD)], [y_fit(1) y_fit(length(y_fit))],'Color',[1 redS(pi2) redS(pi2)])
            
            hidallSTATS.PD{pi2}.brN
            
            tmpMin2 = min(ga2);
            tmpMax2 = max(ga2);
            
            if tmpMin2 < tmpMin
                tmpMin = tmpMin2;
                yMin = tmpMin2;
            else
                yMin = tmpMin;
            end
            
            if tmpMax2 > tmpMax
                tmpMax = tmpMax2;
                yMax = tmpMax2;
            else
                yMax = tmpMax;
            end
            
            ylim([yMin yMax])
            [~ , pval2] = corr(stD,ga2);
            
        end
        
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        ylim([min(yVALs) max(yVALs)])
        
        xticks(linspace(min(stD),max(stD),5))
        xlim([min(stD) max(stD)])

        ylabel('Absolute percent change')
        xlabel('Duration of stimulation (days)')
        
    case 'Thal-Caud' %%% DONE
        
        [~ , thdallSTATS] = subXFStabSEG_v4_HY(1);
        close all
        
        blues = linspace(0,0.9,length(thdallSTATS.PD));
        
        for pi = 1:length(thdallSTATS.PD)
            figure
            
            tmpData = thdallSTATS.PD{pi}.data;
            tmpID = thdallSTATS.PD{pi}.group;
            
            [pval,~,~ ] = ranksum(tmpData(tmpID == 1),tmpData(tmpID == 2));
            
            boxplot(tmpData,tmpID,'Colors',[0 0 0])
            
            yVALs = get(gca,'YTick');
            
            yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
            
            set(gca,'XTick',1:2)
            set(gca,'XTickLabel',{'NoStim','Stim'})
            title([thdallSTATS.PD{pi}.brN ,'  ', num2str(pval)])
            
            scatDat = [[ thdallSTATS.PD{1}.data(thdallSTATS.PD{1}.group == 1) ;...
                thdallSTATS.PD{1}.data(thdallSTATS.PD{1}.group == 2) ],...
                [ones(size(thdallSTATS.PD{1}.data(thdallSTATS.PD{1}.group == 1))) ;...
                repmat(2,size(thdallSTATS.PD{1}.data(thdallSTATS.PD{1}.group == 2)))]];
            
            scatCols = [repmat([0 0 0],size(thdallSTATS.PD{1}.data(thdallSTATS.PD{1}.group == 1))) ;...
                repmat([blues(pi) blues(pi) 1],size(thdallSTATS.PD{1}.data(thdallSTATS.PD{1}.group == 2)))];
            
            hold on
            
            scatter(scatDat(:,2),scatDat(:,1),40,scatCols,'filled')
            
            %% Add p-value
            text(2,max(yVALs) - abs(min(yVALs)/2), sprintf('p = %0.3f',pval))
            
            axis square
            
        end
        close all
        figure
        tmpMin = 0;
        tmpMax = 5;
        for pi2 = 1:length(thdallSTATS.PD)
            
            
            
            %         g1 = hidallSTATS.PD{2}.data(hidallSTATS.PD{2}.group == 1);
            g2 = thdallSTATS.PD{pi2}.data(thdallSTATS.PD{pi2}.group == 2);
            %         ga1 = abs(g1);
            ga2 = abs(g2);
            t = readtable('ob_subj_data.csv');
            tPD = t(ismember(t.cond,'PD'),:);
            tPDg2 = tPD(tPD.groupN == 2,:);
            %         tPDg1 = tPD(tPD.groupN == 1,:);
            
            %             tPDg3 = tPDg2(tPDg2.f_surg_n == gr2IDS,:);
            stD = tPDg2.stimdur;
            %         nstD = tPDg1.stimdur;
            title([thdallSTATS.PD{pi2}.brN ,'  ', num2str(pval)])
            %             plot(stD,ga2,'o','Color',[blues(pi2) blues(pi2) 1])
            hold on
            %             line([min(get(gca,'XTick')) , max(get(gca,'XTick'))],[0 0],'LineStyle','--','Color','k')
            
            coef_fit = polyfit(stD,ga2,1);
            y_fit = polyval(coef_fit,stD);
            
            line([min(stD) max(stD)], [y_fit(1) y_fit(length(y_fit))],'Color',[blues(pi2) blues(pi2) 1])
            
            
            tmpMin2 = min(ga2);
            tmpMax2 = max(ga2);
            
            if tmpMin2 < tmpMin
                yMin = tmpMin2;
            else
                yMin = tmpMin;
            end
            
            if tmpMax2 > tmpMax
                yMax = tmpMax2;
            else
                yMax = tmpMax;
            end
            
            ylim([yMin yMax])
            
            [~ , pval2] = corr(stD,ga2)
            
            
        end
        
    case 'OBFc' %%% DONE
        
        [~ , loballSTATS] = subXFStabHemi_v4_HY(2);
        
        greens = linspace(0,0.7,length(loballSTATS.PD));
        
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
            
            yVALs = get(gca,'YTick');
            
            yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
            
            scatDat = [[ loballSTATS.PD{1}.data(loballSTATS.PD{1}.group == 1) ;...
                loballSTATS.PD{1}.data(loballSTATS.PD{1}.group == 2) ],...
                [ones(size(loballSTATS.PD{1}.data(loballSTATS.PD{1}.group == 1))) ;...
                repmat(2,size(loballSTATS.PD{1}.data(loballSTATS.PD{1}.group == 2)))]];
            
            scatCols = [repmat([0 0 0],size(loballSTATS.PD{1}.data(loballSTATS.PD{1}.group == 1))) ;...
                repmat([greens(pi) 1 greens(pi)],size(loballSTATS.PD{1}.data(loballSTATS.PD{1}.group == 2)))];
            
            hold on
            
            scatter(scatDat(:,2),scatDat(:,1),40,scatCols,'filled')
            
            %% Add p-value
            text(2,max(yVALs) - abs(min(yVALs)/2), sprintf('p = %0.3f',pval))
            
            axis square
            
        end
        close all
        figure
        tmpMin = 0;
        tmpMax = 5;
        for pi2 = 1:length(loballSTATS.PD)
            
            
            
            %         g1 = hidallSTATS.PD{2}.data(hidallSTATS.PD{2}.group == 1);
            g2 = loballSTATS.PD{pi2}.data(loballSTATS.PD{pi2}.group == 2);
            %         ga1 = abs(g1);
            ga2 = abs(g2);
            t = readtable('ob_subj_data.csv');
            tPD = t(ismember(t.cond,'PD'),:);
            tPDg2 = tPD(tPD.groupN == 2,:);
            %         tPDg1 = tPD(tPD.groupN == 1,:);
            
            %             tPDg3 = tPDg2(tPDg2.f_surg_n == gr2IDS,:);
            stD = tPDg2.stimdur;
            %         nstD = tPDg1.stimdur;
            title([loballSTATS.PD{pi2}.brN ,'  ', num2str(pval)])
            %             plot(stD,ga2,'o','Color',[blues(pi2) blues(pi2) 1])
            hold on
            %             line([min(get(gca,'XTick')) , max(get(gca,'XTick'))],[0 0],'LineStyle','--','Color','k')
            
            coef_fit = polyfit(stD,ga2,1);
            y_fit = polyval(coef_fit,stD);
            
            line([min(stD) max(stD)], [y_fit(1) y_fit(length(y_fit))],'Color',[greens(pi2) 1 greens(pi2)])
            
            
            tmpMin2 = min(ga2);
            tmpMax2 = max(ga2);
            
            if tmpMin2 < tmpMin
                yMin = tmpMin2;
            else
                yMin = tmpMin;
            end
            
            if tmpMax2 > tmpMax
                yMax = tmpMax2;
            else
                yMax = tmpMax;
            end
            
            ylim([yMin yMax])
            
            [~ , pval2] = corr(stD,ga2)
            
            
        end
        
        
        
    case 'OB' %%% DONE
        
        %         [~ , OBallSTATS] = subOBVOL_v2_HY(1);
        
        [~ , OBallSTATS] = subOBVOL_v2c_HY(1,1);
        
        [outSTAT2pd,~,~] = ranksum(OBallSTATS.PD{1}.data(OBallSTATS.PD{1}.group == 1) , OBallSTATS.PD{1}.data(OBallSTATS.PD{1}.group == 2));
        
        close all
        boxplot(OBallSTATS.PD{2}.data , OBallSTATS.PD{2}.group)
        
        set(gca,'XTick',1:2)
        set(gca,'XTickLabel',{'NoStim','Stim'})
        
        
        %         set(gca,'YLim',[-0.1 0.1])
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
        [rho1 , pval1] = corr(stD,g2)
        
        figure;
        plot(stD,ga2,'ro')
        hold on
        line([min(get(gca,'XTick')) , max(get(gca,'XTick'))],[0 0],'LineStyle','--','Color','k')
        lsline
        plot(nstD,ga1,'bo')
        [rho2 , pval2] = corr(stD,ga2)
        
        
        
        
        
        
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

