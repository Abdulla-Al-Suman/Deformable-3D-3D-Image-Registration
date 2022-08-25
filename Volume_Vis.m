clc
clear
load('Patient-17CI128128128-MRI.mat')
[d1,d2,d3]=size(V);
[xOC,yOC,zOC] = meshgrid((0:d2-1)-d2/2,(0:d1-1)-d1/2,(0:d3-1)-d3/2);

Parameters5=[0.5000 0.0000 1.0000 0.0000 1.0000 0.0000];
Warped_voleme5 = warp_WindowZchange_3D(V,Parameters5,xOC,yOC,zOC);

Parameters15=[1.5000 0.0000 1.0000 0.0000 1.0000 0.0000];
Warped_voleme15 = warp_WindowZchange_3D(V,Parameters15,xOC,yOC,zOC);

for i=1:3:d3             
    imagesc([V(:,:,i) Warped_voleme5(:,:,i) Warped_voleme15(:,:,i)]), colormap(gray(256));
    axis off
    hold on
    brighten(.2)
    pause
end


figure(1)
vol3d('cdata',V);
view(180,90);
axis equal
axis vis3d

% figure(2)
% vol3d('cdata',Scalled17);
% view(3);
% axis equal
% axis vis3d
% title('0.5')
% 
% figure(3)
% vol3d('cdata',Scalled17s15);
% view(3);
% axis equal
% axis vis3d
% title('1.5')

figure(4)
vol3d('cdata',Warped_voleme5);
view(180,90);
axis equal
axis vis3d
title('0.5 center origin')

figure(5)
vol3d('cdata',Warped_voleme15);
view(180,90);
axis equal
axis vis3d
title('1.5 center origin')
