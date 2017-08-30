function [] = HY_Freesurf_Figure1_map(brArea,cols,alimV)

% HY_Freesurf_multi_MaskSingleTest(13,[1 0 0],0.5)
%%%% NEED THE NUMBERS FOR:
% THALAMUS, CAUDATE, HIPPOCAMPUS
% 49 R Thalamus, 10 L Thalamus

% Hippocampus

% Hippocampus
freeSurf_GM = load_nii('Z:\Yilma_Project\Case_Data\Case_223\Freesurfer\mri\nii\aparc+aseg.nii');
freeSurf_HP = load_nii('Z:\Yilma_Project\Case_Data\Case_223\Freesurfer\mri\nii\lh.hippoSfLabels-T1.v10.FSvoxelSpace.nii');
freeSurf_WM = load_nii('Z:\Yilma_Project\Case_Data\Case_223\Freesurfer\mri\nii\wmparc.nii');
freeSurf_OB = load_nii('Z:\Yilma_Project\Case_Data\Case_223\Freesurfer\mri\nii\wmparc.nii');

hippDAT = {'HATA','parasubiculum','fimbria','hippocampal-fissure','CA3','CA4',...
           'GC-ML-DG','presubiculum','subiculum','Hippocampal_tail',...
           'molecular_layer_HP','CA1'};

hippLabs = [211, 203, 212, 215, 208, 209, 210, 204, 205, 214, 226, 206];

hippInds = ismember(hippDAT,{'CA1','CA3','CA4'});

freeGM = freeSurf_GM.img;
freeHP = freeSurf_HP.img;
freeWM = freeSurf_WM.img;


%  load mydata.mat
%  [node,elem,face]=v2m(mydata,0.5,5,100);
%  plotmesh(node,face)

%% Hippo test



allHP = hippLabs(hippInds);


for i = 1:length(allHP)
    
    colsHP = [rand(1) rand(1) rand(1)];
    tmpHP = allHP(i);
    
    [ gmDAT ] = extract3dobject_FS_v001(freeHP,tmpHP,1,1);
    
    [ bPointsaGM , bboundsaGM ] = FreeSurf_Extract( gmDAT );
    
    hold on
    
    trisurf(bboundsaGM,bPointsaGM(:,2),bPointsaGM(:,1),bPointsaGM(:,3),...
        'Facecolor',colsHP,... % cols
        'Edgecolor','none',...
        'FaceAlpha',0.75)
    
    hold on
    
end


%%

[ gmDAT ] = extract3dobject_FS_v001(freeGM,brArea,1,1);

[ bPointsaGM , bboundsaGM ] = FreeSurf_Extract( gmDAT );


%%

close all

mriLoad = load_nii('Z:\Yilma_Project\Case_Data\Case_223\Freesurfer\mri\nii\brain.nii');


[ brainIm] = Process_MRI(mriLoad);

% eleLoc = 'E:\Dropbox\Publications_Meta\InProgress\MW_STNMap_Methods\Matlab_Code\eleExtract.nii';
% [ rightELE ] = extractDBSLeadPoly2Linear( eleLoc  , 3 , 1.3 , 80, 1);

[brainFlr] = flipBrainML(brainIm);
hold on

% s = slice(1:256, 1:256, 1:256, double(brainFlr), 120, 126, 120);
% shading('interp')
% colormap('gray')
%
% set(s(1), 'alphadata', squeeze(double(brainFlr(:,126,:))), 'facealpha','interp');alim([0 alimV]);
% set(s(2), 'alphadata', squeeze(double(brainFlr(120,:,:))), 'facealpha','interp');alim([0 alimV]);
% set(s(3), 'alphadata', squeeze(double(brainFlr(:,:,120))), 'facealpha','interp');alim([0 alimV]);

s1 = slice(1:256, 1:256, 1:256, double(brainFlr), [], 120, []);
s2 = slice(1:256, 1:256, 1:256, double(brainFlr), 120, [], []);
s3 = slice(1:256, 1:256, 1:256, double(brainFlr), [], [], 100);
shading('interp')
colormap('gray')
%
set(s3, 'alphadata', squeeze(double(brainFlr(:,:,100))), 'facealpha','interp');alim([0 alimV]);
set(s1, 'alphadata', squeeze(double(brainFlr(120,:,:))), 'facealpha','interp');alim([0 alimV]);
set(s2, 'alphadata', squeeze(double(brainFlr(:,120,:))), 'facealpha','interp');alim([0 alimV]);


hold on

trisurf(bboundsaGM,bPointsaGM(:,2),bPointsaGM(:,1),bPointsaGM(:,3),...
    'Facecolor',cols,... % cols
    'Edgecolor','none',...
    'FaceAlpha',0.75)


% [~] = plot3DboundaryElec( rightELE );

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
