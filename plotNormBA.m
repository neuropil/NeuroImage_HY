function [] = plotNormBA()
close all
load('TotalHippoTable.mat')

uniHs = unique(hippoSegTable.HippoArea);
numS = 1:length(uniHs);
for ui = 1:length(uniHs)
    
    tUhi = uniHs{ui};
    
    th = hippoSegTable.Volmm3(contains(hippoSegTable.HippoArea,tUhi));
    thN = th/max(th);
    hold on
    scatter(repmat(numS(ui),size(thN)),thN,10,[rand(1) rand(1) rand(1)])
    
    
end

xticks(numS)
xticklabels(uniHs);
xtickangle(45)
ylabel('Normallized volume')
ylim([0 1])

end