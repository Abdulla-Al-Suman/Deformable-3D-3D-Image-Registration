function [Registered_Moving_Volume,Registered_Parameters] = Registration_Affine_EPD(Moving_Volume,Total_Iteration,LT,UT,d1,d2,d3,x,y,z,edge_list_R,n) 
% clc;
% close all
% clear;
% warning off all;
global m sm mi
%TI=150; %TI=total iteration
m=[1.0000  0.0000  0.0000  0.0000   0.0000  1.0000   0.0000  0.0000  0.0000  0.0000  1.0000  0.0000];%initial parameters
sm=zeros(Total_Iteration,12); %sm= save m
%LT=0.1;UT=0.9;
%tic

%**************loading reference and moving image volume ********
%load('AffineSCV-DS128-V37VV3Manual.mat')
% load('Patient-3CI128128128-MRI.mat')
% R0=V;
% clear V ct_info;
% RR=imresize(CT,.25);
% RR=RR(:,:,11:28);   
% [d11,d22,d33]=size(RR);
% [xs,ys,zs] = meshgrid((0:d22-1)-d22/2,(0:d11-1)-d11/2,(0:d33-1)-d33/2);
% [xn,yn,zn] = meshgrid(-(d22/2-d22/256):(d22/2-d22/256-1+d22/2-d22/256)/256:(d22/2-d22/256-1)-(d22/2-d22/256-1+d22/2-d22/256)/256,-(d11/2-d11/256):(d11/2-d11/256-1+d11/2-d11/256)/256:(d11/2-d11/256-1)-(d11/2-d11/256-1+d11/2-d11/256)/256,-(d33/2-d33/256):(d33/2-d33/256-1+d33/2-d33/256)/256:(d33/2-d33/256-1)-(d33/2-d33/256-1+d33/2-d33/256)/256); %[0:45/256:45-45/256]-45/2
% R0=interp3(xs,ys,zs,RR,xn,yn,zn);
% clear CT ct_info RR xs ys zs xn yn zn d11 d22 d33;

% load('Patient-37CI128128128-MRI.mat')
% I0=V;
% clear V ct_info;
%clear V VV;
% II=imresize(CT,.25);
% II=II(:,:,10:33);
% [d111,d222,d333]=size(II);
% [xs,ys,zs] = meshgrid((0:d222-1)-d222/2,(0:d111-1)-d111/2,(0:d333-1)-d333/2);
% [xn,yn,zn] = meshgrid(-(d222/2-d222/256):(d222/2-d222/256-1+d222/2-d222/256)/256:(d222/2-d222/256-1)-(d222/2-d222/256-1+d222/2-d222/256)/256,-(d111/2-d111/256):(d111/2-d111/256-1+d111/2-d111/256)/256:(d111/2-d111/256-1)-(d111/2-d111/256-1+d111/2-d111/256)/256,-(d333/2-d333/256):(d333/2-d333/256-1+d333/2-d333/256)/256:(d333/2-d333/256-1)-(d333/2-d333/256-1+d333/2-d333/256)/256); %[0:45/256:45-45/256]-45/2
% I0=interp3(xs,ys,zs,II,xn,yn,zn);

%VV=R0; V=I0;
%save Patient_3_37_Cropped_original_MRI-VP37VVP3 V VV
%clear CT ct_info II xs ys zs xn yn zn d111 d222 d333;
%clear V VV D Time;

%%******************** end of data load********************************



% creating edge image of reference volume
%[d1,d2,d3]=size(Fixed_Volume);

% ER = zeros(d1,d2,d3);
% for i=1:d3
% ER(:,:,i) = edge(Fixed_Volume(:,:,i),'canny',[LT UT],1.5);
% end
%[r,c,v] = ind2sub(size(ER),find(ER == 1)); % convert binary edge map into list of pixel positions
%edge_list_R = [c r v];
%DR = ChamferDis3dinf0S(uint32(ER),100);
%clear r c v ER i;
%[n,~] = size(edge_list_R);
%[x,y,z] = meshgrid((0:d2-1)-d2/2,(0:d1-1)-d1/2,(0:d3-1)-d3/2);
%[x,y,z] = meshgrid(1:d2,1:d1,1:d3);

% initialize EPD similarity measure
Smin = inf;
mi=0; %minimum iteration
m0=m;   % m0= if Smin becomes minimum in first iteration
I = warp_AffineEPD_3D(Moving_Volume,m0,x,y,z); %getting colser to reference
%sepd=zeros(1,TI);
for iter = 1:Total_Iteration
    %iter
    [Smin]= update_m_AffineEPD(I,Smin,edge_list_R,iter,d1,d2,d3,n,x,y,z,LT,UT);
    %Smin
    I = warp_AffineEPD_3D(Moving_Volume,m,x,y,z); % place I0 at updated position
end   
if mi==0
    m=m0; %Find m for which epd minimum 
else
    m=sm(mi,:);
end
%mi
Registered_Parameters=m;
Registered_Moving_Volume = warp_AffineEPD_3D(Moving_Volume,m,x,y,z);
clear sm m mi Smin m0 I;
%figure(5), imagesc([sum(R0,3)   sum(I,3)]), colormap(gray) 
%V=I;VV=R0;
%Time=toc;
%save AffineEPD-DS128-V37VV3 V VV m m0 Smin sm sepd TI mi Time 
% figure(1)
% plot(1:TI,sm(:,1),'b',1:TI,sm(:,2),'r',1:TI,sm(:,3),'c',1:TI,sm(:,4),'y',1:TI,sm(:,5),'m',1:TI,sm(:,6),'*',1:TI,sm(:,7),'+',1:TI,sm(:,8),'d',1:TI,sm(:,9),'s',1:TI,sm(:,10),'o',1:TI,sm(:,11),'v',1:TI,sm(:,12),'>')
% figure(2)
% plot(1:TI,sepd(1,:),'g')



