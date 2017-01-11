%% Step 1

fprintf('Running Step 1 \n')

FreeSurf_stats2txt_HY()

fprintf('Step 1 Complete \n')

%% Step 2

fprintf('Running Step 2 \n')

FreeSurf_WMCombine_HY()

fprintf('Step 2 Complete \n')

%% Step 3

fprintf('Running Step 3 \n')

CAT12_Extract_HY()

fprintf('Step 3 Complete \n')

%% Step 4 

fprintf('Running Step 4 \n')

% [allCaseD , allSTATS] = subXcatTab_v1_HY('CTthick');

[allCaseD , allSTATS] = subXcatTab_v2_HY('CTthick');

kruskalwallis(allSTATS.PD{1,1}.data , allSTATS.PD{1,1}.group)
close all
boxplot(allSTATS.PD{1,1}.data , allSTATS.PD{1,1}.group)

pdS1i = allSTATS.PD{1,1}.group == 1;
pdS1 = allSTATS.PD{1,1}.data(pdS1i);
etS1 = allSTATS.ET{1,1}.data;

condGs = [repmat({'PD'},size(pdS1)) ; repmat({'ET'},size(etS1))];
dataGs = [pdS1 ; etS1];
kruskalwallis(dataGs , condGs);
boxplot(dataGs , condGs);

fprintf('Step 4 Complete \n')

%% Step 5

fprintf('Running Step 5 \n')

[allCaseD , allSTATS] = subXFStabCortM_v2_HY();


%% Step 6

[allCaseD , allSTATS] = subXFStabHemi_v2_HY('latOB');




%% Step 7

[allCaseDthal , allSTATSthal] = subXFStabSEG_v2_HY('Thal');

[allCaseDcaud , allSTATScaud] = subXFStabSEG_v2_HY('Caud');

[allCaseDhippo , allSTATShippo] = subXFStabSEG_v2_HY('Hippo');


