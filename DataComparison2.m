
clc
clear 
close all

load('Patient-61CI128128128-MRI.mat')
x1= V;[M,N,Z]=size(x1);
%[xOC,yOC,zOC] = meshgrid((0:N-1)-N/2,(0:M-1)-M/2,(0:Z-1)-Z/2);
load('Patient-3CI128128128-MRI.mat')
x2= V;
load('Patient-45CI128128128-MRI.mat')
x3= V;
% Parameters5=[1.0000 0.0000 1.0000 0.0000 1.5 0.0000];
% x2 = warp_WindowZchange_3D(V,Parameters5,xOC,yOC,zOC);
%x2= V;

% x1 = x1 - min(x1(:)); 
% x2 = x2 - min(x2(:)); 
% x1 = round(255.4*x1/(max(x1(:)))); 
% x2 = round(255.4*x2/(max(x1(:))));
% [~,~,MI] = ent(x1,x1);
% MI
% [MI2,NMI]=MI_GG(x1,x1);
% MI2
% NMI

%load('Patient 25 original MRI.mat');x2= imresize(CT,0.25);
r=0; c=0;EI=0;
for i=1:3:Z             
    EI = edge(x1(:,:,i),'canny',[0.001 0.1],1.5);
    [r,c] = ind2sub(size(EI),find(EI == 1));
    EI = edge(x2(:,:,i),'canny',[0.001 0.1],1.5);
    [rr,cc] = ind2sub(size(EI),find(EI == 1));
    EI = edge(x3(:,:,i),'canny',[0.04 0.2],1.5);
    [rrr,ccc] = ind2sub(size(EI),find(EI == 1));
    imagesc([x1(:,:,i) x2(:,:,i) x3(:,:,i)]), colormap(gray(256));
    %imagesc(x2(:,:,i)), colormap(gray(256));
    axis off
    hold on
    scatter(c,r,1,'r','LineWidth',2)
    scatter(cc+128,rr,1,'r','LineWidth',2)
    scatter(ccc+256,rrr,1,'r','LineWidth',2)
    brighten(.2)
    pause
    r=0; c=0;EI=0;
end



