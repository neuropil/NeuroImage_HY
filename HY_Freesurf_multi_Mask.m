%%%

freeSurfdat1 = load_nii('E:\Dropbox\PowerPoint_Meta\MHA_talk2016\ExtraImages\mriNII\c202\aparc+aseg.nii');
freeSurfdat2 = load_nii('E:\Dropbox\PowerPoint_Meta\MHA_talk2016\ExtraImages\mriNII\c204\aparc+aseg.nii');
freeSurfdat3 = load_nii('E:\Dropbox\PowerPoint_Meta\MHA_talk2016\ExtraImages\mriNII\c205\aparc+aseg.nii');
freeSurfdat4 = load_nii('E:\Dropbox\PowerPoint_Meta\MHA_talk2016\ExtraImages\mriNII\c212\aparc+aseg.nii');
freeSurfdat5 = load_nii('E:\Dropbox\PowerPoint_Meta\MHA_talk2016\ExtraImages\mriNII\c213\aparc+aseg.nii');
freeSurfdat6 = load_nii('E:\Dropbox\PowerPoint_Meta\MHA_talk2016\ExtraImages\mriNII\c215\aparc+aseg.nii');
freeSurfdat7 = load_nii('E:\Dropbox\PowerPoint_Meta\MHA_talk2016\ExtraImages\mriNII\c216\aparc+aseg.nii');
freeSurfdat8 = load_nii('E:\Dropbox\PowerPoint_Meta\MHA_talk2016\ExtraImages\mriNII\c218\aparc+aseg.nii');
freeSurfdat9 = load_nii('E:\Dropbox\PowerPoint_Meta\MHA_talk2016\ExtraImages\mriNII\c220\aparc+aseg.nii');
freeSurfdat10 = load_nii('E:\Dropbox\PowerPoint_Meta\MHA_talk2016\ExtraImages\mriNII\c223\aparc+aseg.nii');

freeSim1 = freeSurfdat1.img;
freeSim2 = freeSurfdat2.img;
freeSim3 = freeSurfdat3.img;
freeSim4 = freeSurfdat4.img;
freeSim5 = freeSurfdat5.img;
freeSim6 = freeSurfdat6.img;
freeSim7 = freeSurfdat7.img;
freeSim8 = freeSurfdat8.img;
freeSim9 = freeSurfdat9.img;
freeSim10 = freeSurfdat10.img;

allFREES = {freeSim1 , freeSim2 , freeSim3 , freeSim4 , freeSim5 ,...
    freeSim6 , freeSim7 , freeSim8 , freeSim9 , freeSim10};


%%

% randSelBA = randperm(100,25);
bPointsa = cell(length(allFREES),1);
bboundsa = cell(length(allFREES),1);
maxZ = zeros(length(allFREES),1);
minY = zeros(length(allFREES),1);
for i = 1:length(allFREES)
    
    
    tFree = allFREES{i};
    
    [ freeDat ] = extract3dobject_FS_v001(tFree,49,1,0);
    
    [ bPointsa{i} , bboundsa{i} ] = FreeSurf_Extract( freeDat );
    
    maxZ(i) = max(bPointsa{i}(:,3));
    minY(i) = min(bPointsa{i}(:,2));
    
end

maxZ2 = max(maxZ) - maxZ;
minY2 = minY - min(minY);



%%



figure
blobP = 0;
for i = 1:length(allFREES)
    
    
    tFree = allFREES{i};
    
    if i == 1
        blobP = 0;
    else
        blobP = blobP + 20;
    end
    
    cols = [1 0.2 0.1];
    
    [ freeDat ] = extract3dobject_FS_v001(tFree,49,1,0);
    
    [ blobPoints , blobBounds ] = FreeSurf_Extract( freeDat );
    blobPoints(:,1) = blobPoints(:,1) + blobP;
    blobPoints(:,3) = blobPoints(:,3) + maxZ2(i);
    blobPoints(:,2) = blobPoints(:,2) - minY2(i);
    hold on
    
    trisurf(blobBounds,blobPoints(:,2),blobPoints(:,1),blobPoints(:,3),...
        'Facecolor',cols,...
        'Edgecolor','none',...
        'FaceAlpha',0.5)
    
    
    xlim([80 170]);
    ylim([100 320]);
    zlim([50 200]);
    
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        set(gca,'ZTickLabel',[])
    
    
    %     pause
end


%%

close all

mriLoad = load_nii('E:\Dropbox\PowerPoint_Meta\MHA_talk2016\ExtraImages\mriNII\c205\brain.nii');


[ brainIm] = Process_MRI(mriLoad);


% [brainFlr] = flipBrainML(brainIm);
[brainFlr] = flipBrainML(brainIm);
figure;

s = slice(1:256, 1:256, 1:256, double(brainFlr), 120, 120, 115);
shading('interp')
colormap('gray')

set(s(1), 'alphadata', squeeze(double(brainFlr(:,120,:))), 'facealpha','interp');alim([0 0.2]);
set(s(2), 'alphadata', squeeze(double(brainFlr(120,:,:))), 'facealpha','interp');alim([0 0.2]);
set(s(3), 'alphadata', squeeze(double(brainFlr(:,:,115))), 'facealpha','interp');alim([0 0.8]);

hold on

trisurf(bboundsa{3},bPointsa{3}(:,2),bPointsa{3}(:,1),bPointsa{3}(:,3),...
    'Facecolor',cols,...
    'Edgecolor','none',...
    'FaceAlpha',0.5)

xlim([20 200]);
ylim([50 200]);
zlim([20 200]);

set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
set(gca,'ZTickLabel',[])
