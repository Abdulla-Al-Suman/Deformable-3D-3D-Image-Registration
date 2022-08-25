function Vol = ptCloud2vol(X,Y,Z,muscle_number)


V = zeros(128,128,128);

for i = 1:length(X{muscle_number})
    V(round(X{muscle_number}(i)),round(Y{muscle_number}(i)),round(Z{muscle_number}(i))) = 1;
end

% figure(1);
% vol3d('cdata',V,'texture','3D');
% axis vis3d
% axis equal
% view(0,0)

FV = isosurface(V,0.2);

% figure(2); patch(FV,'FaceColor',       [0.8 0.8 1.0], ...
%          'EdgeColor',       'none',        ...
%          'FaceLighting',    'gouraud',     ...
%          'AmbientStrength', 0.15);
% % Add a camera light, and tone down the specular highlighting
% camlight('headlight');
% material('dull');
% view(0,0)

Vf = voxelise(round(max(FV.vertices(:,1)))-round(min(FV.vertices(:,1)))+1,round(max(FV.vertices(:,2)))-round(min(FV.vertices(:,2)))+1,round(max(FV.vertices(:,3)))-round(min(FV.vertices(:,3)))+1,FV,'xyz');

% figure;
% vol3d('cdata',Vf,'texture','3D');
% axis vis3d
% axis equal
% view(0,90)

Vol = zeros(128,128,128);
Vol(round(min(FV.vertices(:,1))):round(max(FV.vertices(:,1))),round(min(FV.vertices(:,2))):round(max(FV.vertices(:,2))),round(min(FV.vertices(:,3))):round(max(FV.vertices(:,3)))) = Vf;

for k = 1:128
    Vol(:,:,k) = imfill(imclose(squeeze(Vol(:,:,k)),strel('disk',2)),'holes');
%     figure(5); imagesc(Vol(:,:,k))
%     pause
end

