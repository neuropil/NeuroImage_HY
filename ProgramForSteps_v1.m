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

[allCaseD , allSTATS] = subXcatTab_v1_HY('CTthick');

kruskalwallis(allSTATS{1,1}.data , allSTATS{1,1}.group)

boxplot(allSTATS{1,1}.data , allSTATS{1,1}.group)

fprintf('Step 4 Complete \n')