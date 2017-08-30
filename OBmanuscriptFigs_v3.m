

% This function


function [outSTAT1pd, outSTAT2pd, outSTAT1etpd] = OBmanuscriptFigs_v3(ID)

switch ID
    
    case 'Hippocampus' %%% DONE
        
        [~ , hidallSTATS, ids] = subXhippo_HYb(2,1,1);
        gr2IDS = transpose([ids{2,:}]);
        redS = linspace(0,0.7,length(hidallSTATS.PD));
        
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
            
            cd('C:\Users\johna\Dropbox\Publications_Meta\InProgress\Yilma_OB_DBS\images\NewFigs08032017')
            
            fname = ['Figure1_',hidallSTATS.PD{pi}.brN,'.pdf'];
            
            saveas(gcf, fname)
            
        end
        close all
        figure
        
        tmpMin = 0;
        tmpMax = 5;
        pVals = zeros(1,length(hidallSTATS.PD));
        rhos = zeros(1,length(hidallSTATS.PD));
        for pi2 = 1:length(hidallSTATS.PD)
            
            cd('C:\Users\johna\Dropbox\Publications_Meta\InProgress\Yilma_OB_DBS\DATA')
            
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
            
            [rhos(pi2) , pVals(pi2)] = corr(stD,ga2);
            
            if rhos(pi2) < 0
                
                line([min(stD) max(stD)], [y_fit(1) y_fit(length(y_fit))],'Color',[1 redS(pi2) redS(pi2)])
                
            else
                
                line([min(stD) max(stD)], [y_fit(length(y_fit)) y_fit(1)],'Color',[1 redS(pi2) redS(pi2)])
                
            end
            
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
            
            
        end
        
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        ylim([min(yVALs) max(yVALs)])
        
        xticks(linspace(min(stD),max(stD),5))
        xlim([min(stD) max(stD)])
        
        ylabel('Absolute percent change')
        xlabel('Duration of stimulation (days)')
        
        yT = max(yVALs);
        xT = get(gca,'XTick');
        xT2 = xT(3);
        for ti = 1:length(hidallSTATS.PD)
            
            hidallSTATS.PD{ti}.brN = replace(hidallSTATS.PD{ti}.brN,'_',' ');
            
            text(xT2,yT - 1,hidallSTATS.PD{ti}.brN,'Color',[1 redS(ti) redS(ti)])
            text(xT2 + 200,yT - 1,sprintf('%0.2f',pVals(ti)),'Color','k')
            text(xT2 + 250,yT - 1,sprintf('%0.2f',rhos(ti)),'Color','k')
            
            yT = yT - 1;
        end
        
        cd('C:\Users\johna\Dropbox\Publications_Meta\InProgress\Yilma_OB_DBS\images\NewFigs08032017')
        
        
        axis square
        fname = 'Figure1_HippoCorr.pdf';
        
        saveas(gcf, fname)
        
        close all
        
        
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
            
            cd('C:\Users\johna\Dropbox\Publications_Meta\InProgress\Yilma_OB_DBS\images\NewFigs08032017')
            
            fname = ['Figure2_',loballSTATS.PD{pi}.brN,'.pdf'];
            
            saveas(gcf, fname)
            
        end
        close all
        figure
        tmpMin = 0;
        tmpMax = 5;
        pVals = zeros(1,length(loballSTATS.PD));
        rhos = zeros(1,length(loballSTATS.PD));
        for pi2 = 1:length(loballSTATS.PD)
            
            cd('C:\Users\johna\Dropbox\Publications_Meta\InProgress\Yilma_OB_DBS\DATA')
            
            g2 = loballSTATS.PD{pi2}.data(loballSTATS.PD{pi2}.group == 2);
            
            ga2 = abs(g2);
            t = readtable('ob_subj_data.csv');
            tPD = t(ismember(t.cond,'PD'),:);
            tPDg2 = tPD(tPD.groupN == 2,:);
            stD = tPDg2.stimdur;
            
            title('Olfactory Cortex')
            
            hold on
            
            coef_fit = polyfit(stD,ga2,1);
            y_fit = polyval(coef_fit,stD);
            [rhos(pi2) , pVals(pi2)] = corr(stD,ga2);
            if rhos(pi2) < 0
                
                line([min(stD) max(stD)], [y_fit(1) y_fit(length(y_fit))],'Color',[greens(pi2) 1 greens(pi2)])
                
            else
                
                line([min(stD) max(stD)], [y_fit(length(y_fit)) y_fit(1)],'Color',[greens(pi2) 1 greens(pi2)])
                
            end
            
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
            
            
        end
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        ylim([min(yVALs) max(yVALs)])
        
        xticks(linspace(min(stD),max(stD),5))
        xlim([min(stD) max(stD)])
        
        ylabel('Absolute percent change')
        xlabel('Duration of stimulation (days)')
        
        yT = max(yVALs);
        xT = get(gca,'XTick');
        xT2 = xT(3);
        for ti = 1:length(loballSTATS.PD)
            
            loballSTATS.PD{ti}.brN = replace(loballSTATS.PD{ti}.brN,'_',' ');
            
            text(xT2,yT - 2,loballSTATS.PD{ti}.brN,'Color',[greens(ti) 1 greens(ti)])
            text(xT2 + 200,yT - 2,sprintf('%0.2f',pVals(ti)),'Color','k')
            text(xT2 + 300,yT - 2,sprintf('%0.2f',rhos(ti)),'Color','k')
            
            yT = yT - 2;
        end
        cd('C:\Users\johna\Dropbox\Publications_Meta\InProgress\Yilma_OB_DBS\images\NewFigs08032017')
        axis square
        
        fname = 'Figure2_OC_Corr.pdf';
        
        saveas(gcf, fname)
        
    case 'OB' %%% DONE
        
        %         [~ , OBallSTATS] = subOBVOL_v2_HY(1);
        
        [~ , OBallSTATS] = subOBVOL_v2c_HY(1,1);
        
        [pval,~,~] = ranksum(OBallSTATS.PD{1}.data(OBallSTATS.PD{1}.group == 1) , OBallSTATS.PD{1}.data(OBallSTATS.PD{1}.group == 2));
        
        close all
        boxplot(OBallSTATS.PD{2}.data , OBallSTATS.PD{2}.group)
        
        set(gca,'XTick',1:2)
        set(gca,'XTickLabel',{'NoStim','Stim'})
        
        title(['OB  ', num2str(pval)])
        
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        
        scatDat = [[ OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 1) ;...
            OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 2) ],...
            [ones(size(OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 1))) ;...
            repmat(2,size(OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 2)))]];
        
        scatCols = [repmat([0 0 0],size(OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 1))) ;...
            repmat([0 0 1],size(OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 2)))];
        
        hold on
        
        scatter(scatDat(:,2),scatDat(:,1),40,scatCols,'filled')
        axis square
        
        hostname = char( getHostName( java.net.InetAddress.getLocalHost ) );
        
        if strcmp(hostname,'JAT-PC')
            
            cd('E:\Dropbox\Publications_Meta\InProgress\Yilma_OB_DBS\images\NewFigs08032017')
            
        else
            
            cd('C:\Users\johna\Dropbox\Publications_Meta\InProgress\Yilma_OB_DBS\images\NewFigs08032017')
            
        end
        
        
        
        fname = 'Figure3_OB.pdf';
        
        saveas(gcf, fname)
        
        close all
        
        if strcmp(hostname,'JAT-PC')
            
            cd('E:\Dropbox\Publications_Meta\InProgress\Yilma_OB_DBS\DATA')
            
        else
            
            cd('C:\Users\johna\Dropbox\Publications_Meta\InProgress\Yilma_OB_DBS\DATA')
            
        end
        
        %         g1 = OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 1);
        g2 = OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 2);
        g1 = OBallSTATS.PD{2}.data(OBallSTATS.PD{2}.group == 1);
        %         ga1 = abs(g1);
        ga2 = abs(g2);
        ga1 = abs(g1);
        t = readtable('ob_subj_data.csv');
        obb = readtable('obBehv.csv');
        tPD = t(ismember(t.cond,'PD'),:);
        tPDg2 = tPD(tPD.groupN == 2,:);
        tPDg1 = tPD(tPD.groupN == 1,:);
        %         tPDg1 = tPD(tPD.groupN == 1,:);
        stD = tPDg2.stimdur;
        stD1 = tPDg1.stimdur;
        %         nstD = tPDg1.stimdur;
        
        %         figure;
        %         plot(stD,g2,'ro')
        %         hold on
        %         line([min(get(gca,'XTick')) , max(get(gca,'XTick'))],[0 0],'LineStyle','--','Color','k')
        %         lsline
        %         [rho1 , pval1] = corr(stD,g2)
        
        figure;
        plot(stD,ga2,'bo')
        hold on
        
        coef_fit = polyfit(stD,ga2,1);
        y_fit = polyval(coef_fit,stD);
        
        line([min(stD) max(stD)], [y_fit(length(y_fit)) y_fit(1)],'Color','b')
        
        yVALs = get(gca,'YTick');
        
        yticks([min(yVALs) round((max(yVALs) - min(yVALs))/2,3)+min(yVALs)  max(yVALs)])
        ylim([min(yVALs) max(yVALs)])
        
        xticks(linspace(min(stD),max(stD),5))
        xlim([min(stD) max(stD)])
        
        ylabel('Absolute percent change')
        xlabel('Duration of stimulation (days)')
        
        [rho1 , pval2] = corr(stD,ga2);
        
        axis square
        
        yT = max(yVALs);
        xT = get(gca,'XTick');
        xT2 = xT(3);
        
        
        text(xT2 + 200,yT - 2,sprintf('p = %0.2f',pval2),'Color','k')
        text(xT2 + 200,yT - 5,sprintf('r = %0.2f',rho1),'Color','k')
        
        
        if strcmp(hostname,'JAT-PC')
            
            cd('E:\Dropbox\Publications_Meta\InProgress\Yilma_OB_DBS\images\NewFigs08032017')
            
        else
            
            cd('C:\Users\johna\Dropbox\Publications_Meta\InProgress\Yilma_OB_DBS\images\NewFigs08032017')
            
        end
        
        fname = 'Figure3_OB_Corr.pdf';
        
        axis square
        
        saveas(gcf, fname)
        close all
        
        plot(obb.stimbetw,obb.UPSIT,'ko')
        [rho_vu , pval_vu] = corr(obb.stimbetw,obb.UPSIT);
        
        coef_fit_vu = polyfit(obb.stimbetw,obb.UPSIT,1);
        y_fit_vu = polyval(coef_fit_vu,obb.stimbetw);
        hold on
        line([min(obb.stimbetw) max(obb.stimbetw)], [y_fit_vu(length(y_fit_vu)) y_fit_vu(1)],'Color','k','LineWidth',2.5)
        
        xVals1 = get(gca,'XTick');
        yVals1 = get(gca,'YTick');
        xticks(linspace(min(xVals1),max(xVals1),3))
        yticks(linspace(min(yVals1),max(yVals1),3))
        
        ylabel('Pre/Post stimulation change in UPSIT')
        xlabel('Duration of stimulation before test (days)')
        axis square
        title(sprintf('p = %0.2f r = %0.2f',rho_vu , pval_vu))
        
        obbInd2 = ismember(tPDg2.f_surg_n,obb.first);
        obbTab2 = tPDg2(obbInd2,:);
        obbTab2.OBvol = g2(obbInd2);
        obbBind2 = ismember(obb.first,obbTab2.f_surg_n);
        obbTab2.Upsit = obb.UPSIT(obbBind2);
        obbTab2.stimB = obb.stimbetw(obbBind2);
        
        obbInd1 = ismember(tPDg1.f_surg_n,obb.first);
        obbTab1 = tPDg1(obbInd1,:);
        obbTab1.OBvol = g1(obbInd1);
        obbBind1 = ismember(obb.first,obbTab1.f_surg_n);
        obbTab1.Upsit = obb.UPSIT(obbBind1);
        obbTab1.stimB = obb.stimbetw(obbBind1);
        
        figure;
        %         nstimB = (obbTab.stimB - min(obbTab.stimB)) ./ (max(obbTab.stimB) - min(obbTab.stimB));
        %         nOBvol = (obbTab.OBvol - min(obbTab.OBvol)) ./ (max(obbTab.OBvol) - min(obbTab.OBvol));
%         subplot(1,2,1)
%         plot(obbTab.stimB,obbTab.Upsit,'bo')
%         [rho_su , pval_su] = corr(obbTab.stimB,obbTab.Upsit);
%         
%         coef_fit_su = polyfit(obbTab.stimB,obbTab.Upsit,1);
%         y_fit_su = polyval(coef_fit_su,obbTab.stimB);
%         hold on
%         line([min(obbTab.stimB) max(obbTab.stimB)], [y_fit_su(length(y_fit_su)) y_fit_su(1)],'Color','b','LineWidth',2.5)
%         
%         title(sprintf('r = %0.2f , p = %0.2f',rho_su,pval_su))
        
%         xticks(linspace(min(obbTab.stimB),max(obbTab.stimB),5))
%         xlim([min(obbTab.stimB) max(obbTab.stimB)])
%         
%         ylabel('Change in UPSIT Pre/Post DBS')
%         xlabel('Duration of stimulation')
        
%         subplot(1,2,2)
% 
%         OBdat = [obbTab1.OBvol ; obbTab2.OBvol];
%         
%         OBgr = [repmat(1,size(obbTab1.OBvol)) ; repmat(2,size(obbTab2.OBvol))]
%         
%         boxplot(OBdat,OBgr);
%         [a,b,c] = ranksum(obbTab1.OBvol , obbTab2.OBvol);
        
        figure;

        plot(obbTab2.OBvol,obbTab2.Upsit,'ro')
        [rho_vu , pval_vu] = corr(obbTab2.OBvol,obbTab2.Upsit);
        
        coef_fit_vu = polyfit(obbTab2.OBvol,obbTab2.Upsit,1);
        y_fit_vu = polyval(coef_fit_vu,obbTab2.OBvol);
        hold on
        line([min(obbTab2.OBvol) max(obbTab2.OBvol)], [y_fit_vu(length(y_fit_vu)) y_fit_vu(1)],'Color','r','LineWidth',2.5)
        
        
        plot(obbTab1.OBvol,obbTab1.Upsit,'bo')
        [rho_vu1 , pval_vu1] = corr(obbTab1.OBvol,obbTab1.Upsit);
        
        coef_fit_vu1 = polyfit(obbTab1.OBvol,obbTab1.Upsit,1);
        y_fit_vu1 = polyval(coef_fit_vu1,obbTab1.OBvol);
        hold on
        line([min(obbTab1.OBvol) max(obbTab1.OBvol)], [y_fit_vu1(length(y_fit_vu1)) y_fit_vu1(1)],'Color','b','LineWidth',2.5)
        
        
        
        title(sprintf('r = %0.2f , p = %0.2f',rho_vu,pval_vu))
        
        xticks(round(linspace(3,max(obbTab2.OBvol),5)))
        xlim([3 max(obbTab2.OBvol)])
        
        xlabel('Change in OB volume post stimulation')
        ylabel('Change in UPSIT Pre/Post DBS')
        
        axis square
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

