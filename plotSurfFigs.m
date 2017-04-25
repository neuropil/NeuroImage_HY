
% directory = 'E:\Dropbox\ChungMooMeta';
% SPHARMconstruct(directory,85);
% 
% 
% %%
% 
% [tri,coord,nbr,normal]=mni_getmesh('outersurface.obj');
% trisurf(tri,coord(1,:),coord(2,:),coord(3,:))
% 
% 
% %%
% 
%  [node,elem]=readasc('lh.orig.asc');
% %  read FreeSurfer ASC mesh format
% %  input:
% %       fname: name of the asc file
% %  output:
% %       node: node positions of the mesh
% %       elem: element list of the mesh
% 
% 
% %%
% 
%  saveasc(node,elem,'lhorig.obj')
% %  save a surface mesh to FreeSurfer ASC mesh format
% %  input:
% %       v: input, surface node list, dimension (nn,3)
% %       f: input, surface face element list, dimension (be,3)
% %       fname: output file name
% 
% %%
% % vertface2obj(v,f,'gyroid.obj')
% vertface2obj(node,elem,'gyroid.obj')
% 
% %%
% 
% [tri,coord,nbr,normal]=mni_getmesh('gyroid.obj');
% trisurf(tri,coord(1,:),coord(2,:),coord(3,:))

%%

[vertices, faces] = read_surf('lh.pial');

%%
surf.vertices=vertices;
surf.faces=faces+1;

figure_wire(surf,'white','none');

%% 
% surf.vertices=squeeze(surf.vertices(1,:,:))';
% surf.faces=tri;
surf2 = reducepatch(surf,0.75);
figure; figure_wire(surf2, 'r', 'none');

%%

figure_patch(surf2,'k',0.5);
hold on
figure_wire(surf2, 'w', 'none');


%%


dir='lh.sphere'
[vertices2, faces2] = read_surf(dir);
sphere.vertices=vertices2;
sphere.faces=faces2+1;
figure_wire(sphere,'white','none');

%%

[theta varphi]=EULERangles(sphere);
figure;
figure_trimesh(surf,theta,'rwb');
figure;
figure_trimesh(surf,varphi,'rwb'); 

%%

[surf_smooth, fourier]=SPHARMsmooth2(surf,sphere,40,0);

%%
figure_wire(surf_smooth,'red','none');

%%

[vertices, faces] = read_surf('lh.white');

surf.vertices=vertices;
surf.faces=faces+1;

figure_wire(surf,'white','none');

%%

[theta varphi]=EULERangles(sphere);
figure;
figure_trimesh(surf,theta,'rwb');
figure;
figure_trimesh(surf,varphi,'rwb'); 

%%

[surf_smooth, fourier]=SPHARMsmooth2(surf,sphere,40,0);

%
figure_wire(surf_smooth,'red','none');

%% Case 229

[vertices, faces] = read_surf('lh.pial');
surf.vertices=vertices;
surf.faces=faces+1;

dir='lh.sphere';
[vertices2, faces2] = read_surf(dir);
sphere.vertices=vertices2;
sphere.faces=faces2+1;
[theta varphi]=EULERangles(sphere);
[surf_smooth, fourier]=SPHARMsmooth2(surf,sphere,40,0);
figure_wire(surf_smooth,'m','none');


hold on

[vertices, faces] = read_surf('rh.pial');
surf.vertices=vertices;
surf.faces=faces+1;

dir='rh.sphere';
[vertices2, faces2] = read_surf(dir);
sphere.vertices=vertices2;
sphere.faces=faces2+1;
[theta varphi]=EULERangles(sphere);
[surf_smooth, fourier]=SPHARMsmooth2(surf,sphere,40,0);
% figure;
figure_wire(surf_smooth,[0.5 0.5 0.5],'none');
