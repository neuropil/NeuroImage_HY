function [] = HY_Freesurf_multi_MaskSingleTest(brArea,cols,alimV)

freeSurfdat1 = load_nii('E:\Dropbox\Publications_Meta\InProgress\MW_STNMap_Methods\Matlab_Code\aparc+aseg_E_DTI.nii');

freeSim1 = freeSurfdat1.img;


%%

[ freeDat ] = extract3dobject_FS_v001(freeSim1,brArea,1,0);

[ bPointsa , bboundsa ] = FreeSurf_Extract( freeDat );


%%

close all

mriLoad = load_nii('E:\Dropbox\Publications_Meta\InProgress\MW_STNMap_Methods\Matlab_Code\brain.nii');


[ brainIm] = Process_MRI(mriLoad);

eleLoc = 'E:\Dropbox\Publications_Meta\InProgress\MW_STNMap_Methods\Matlab_Code\eleExtract.nii';
[ rightELE ] = extractDBSLeadPoly2Linear( eleLoc  , 3 , 1.3 , 80, 1);
% [ rightELE ] = extractDBSLeadPolyWK_v1( eleLoc  , 2 , [0.5 0.5 0.5] , 1.3 , 80 , 'right', 1);
% [boundaryObjR,circleContainerR] = plot3DboundaryElec_DK_v3( rightELE );


% [brainFlr] = flipBrainML(brainIm);
[brainFlr] = flipBrainML(brainIm);
hold on

% s = slice(1:256, 1:256, 1:256, double(brainFlr), 120, 126, 120);
% shading('interp')
% colormap('gray')
%
% set(s(1), 'alphadata', squeeze(double(brainFlr(:,126,:))), 'facealpha','interp');alim([0 alimV]);
% set(s(2), 'alphadata', squeeze(double(brainFlr(120,:,:))), 'facealpha','interp');alim([0 alimV]);
% set(s(3), 'alphadata', squeeze(double(brainFlr(:,:,120))), 'facealpha','interp');alim([0 alimV]);

s = slice(1:256, 1:256, 1:256, double(brainFlr), [], 120, []);
shading('interp')
colormap('gray')
%
% set(s, 'alphadata', squeeze(double(brainFlr(:,126,:))), 'facealpha','interp');alim([0 alimV]);
set(s, 'alphadata', squeeze(double(brainFlr(120,:,:))), 'facealpha','interp');alim([0 alimV]);
% set(s(3), 'alphadata', squeeze(double(brainFlr(:,:,120))), 'facealpha','interp');alim([0 alimV]);


hold on

trisurf(bboundsa,bPointsa(:,2),bPointsa(:,1),bPointsa(:,3),...
    'Facecolor',cols,...
    'Edgecolor','none',...
    'FaceAlpha',0.75)


[~] = plot3DboundaryElec( rightELE );

xlim([20 220]);
ylim([50 200]);
zlim([20 200]);

set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
set(gca,'ZTickLabel',[])

set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'ZTick',[])

set(gca,'XColor','none')
set(gca,'YColor','none')
set(gca,'ZColor','none')

set(gca,'Color','none')

set(gca, 'View', [0 0])

axis square
