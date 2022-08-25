function [Registered_Moving_Volume,Registered_Parameters] = Registration_FFD_EPD(Moving_Volume,Total_Iteration,LT,UT,d1,d2,d3,x,y,z,edge_list_R,n,NB,dxdm)
% clc;
% close all
% clear;
% warning off all;
global m sm mi
%tic
%************** START loading reference and moving image volume ********
%load('FFD5432-DS128-V37VV3Manual.mat')
%load('Patient 30 original MRI.mat') %RR=R(:,:,11:34); pt 30 ROI %III=II(:,:,10:34); pt56 ROI
% RR=single(imresize(CT,.25));
% %  RR=VV; %VV(:,:,80:180);
% [d11,d22,d33]=size(RR);
% [xs,ys,zs] = meshgrid((0:d22-1)-d22/2,(0:d11-1)-d11/2,(0:d33-1)-d33/2);
% [xn,yn,zn] = meshgrid(-(d22/2-d22/256):(d22/2-d22/256-1+d22/2-d22/256)/256:(d22/2-d22/256-1)-(d22/2-d22/256-1+d22/2-d22/256)/256,-(d11/2-d11/256):(d11/2-d11/256-1+d11/2-d11/256)/256:(d11/2-d11/256-1)-(d11/2-d11/256-1+d11/2-d11/256)/256,-(d33/2-d33/256):(d33/2-d33/256-1+d33/2-d33/256)/256:(d33/2-d33/256-1)-(d33/2-d33/256-1+d33/2-d33/256)/256); %[0:45/256:45-45/256]-45/2
% R0=interp3(xs,ys,zs,RR,xn,yn,zn);
%R0=VV;
%clear V ct_info;

%load('Patient-37CI128128128-MRI.mat')
% load('Patient 56 original MRI.mat')
% II=single(imresize(CT,.25));
% % II=V;
% I0=interp3(xs,ys,zs,II,xn,yn,zn);
%I0=V;
%clear V VV;
%clear CT ct_info II RR xs ys zs xn yn zn d11 d22 d33; 
%%******************** end of data load********************************

 %************** START Initialization ********
%[d1,d2,d3]=size(R0);
% TI=100; %TI=total iteration
% NSubband = 1;
% LT=0.01; UT=0.1;
% jj=1; % original spline support=4, jj=6 means support=256, 256/4=64=2^6 
% NB=((2^-jj)*d1+1)*NSubband; % NB= NUMBER OF BASIS FUNCTION
%NB=((2^-jj)*256-1)*NSubband; % NB= NUMBER OF BASIS FUNCTION
m=zeros(1,3*NB);
%m=zeros(1,3*NB,'single');
sm=zeros(Total_Iteration,3*NB); %sm= save m, iteration,no. of parameters

%[x,y,z] = meshgrid((0:d2-1)-d2/2,(0:d1-1)-d1/2,(0:d3-1)-d3/2);
%[x,y,z] = meshgrid((0:d2-1),(0:d1-1),(0:d3-1));
%[x,y,z] = meshgrid(1:d2,1:d1,1:d3);
%************** END Initialization ********


%******** START Creating edge image of reference volume **********

% ER = zeros(d1,d2,d3);
% for i=1:d3
%     ER(:,:,i) = edge(R0(:,:,i),'canny',[LT UT],1.5); % need to select threshold carefully
% end
% % DR = ChamferDis3dinfv2(uint32(ER),inf);
% % %DR = imChamferDistance3d(ER);
% % %DRR = ChamferDis3dinf(uint32(ER),inf);
% % figure(1), imagesc(ER(:,:,2)), colormap(gray)
% % figure(2), imagesc(DR(:,:,2)), colormap(gray)
% %figure(3), imagesc(DRR(:,:,2)), colormap(gray)
% 
% [r,c,v] = ind2sub(size(ER),find(ER == 1)); % convert binary edge map into list of pixel positions
% edge_list_R = [c r v];
% [n,~] = size(edge_list_R);
% clear r c v ER i;

%******** END Creating edge image of reference volume **********

%********** Basis Function Generation Start**********
% dxdm = zeros(d1,d2,d3,NB);
% % TensorP2 = zeros(d1,d2,d3,'single');
% % TensorP3 = zeros(d1,d2,d3,'single');
% %X=0:(d2-1);%Y=0:255;Z=0:255;
% X=1:d2;
% ii=0;
% for k=-((2^-jj)*d1)/2:((2^-jj)*d1)/2  % control points, translation index
% %for k=-(((2^-jj)*d1)/2-1):(((2^-jj)*d1)/2-1)
%     TensorP2 = zeros(d1,d2,d3);
%     ii=ii+1;
%     Spline = SplineST(X,jj,k);
%     %SplineWavelet = Spline_wavelet2(X,jj,k);
%     KronProduct1 = kron(Spline,Spline');
%     %     KronProduct2 = kron(SplineWavelet,Spline');
%     for i=1:d3
%        % OrderPairs1=KronProduct1*Spline(i);
%         OrderPairs2=KronProduct1*Spline(i);
% %         OrderPairs3=KronProduct2*Spline(i);
%         %TensorP1(:,:,i)=OrderPairs1;
%         TensorP2(:,:,i)=OrderPairs2;
% %         TensorP3(:,:,i)=OrderPairs3;
% clear OrderPairs2;
%     end
%     %dxdm(:,:,:,ii)=(2^-jj)*SM;
%     dxdm(:,:,:,ii)=2^-(jj)*TensorP2;
% %     dxdm(:,:,:,ii+((2^-jj)*256+1))=2^-jj*TensorP2;
% %     dxdm(:,:,:,ii+((2^-jj)*256+1)*(NSubband-1))=2^-jj*TensorP3;
%     clear KronProduct1 i TensorP2 Spline;
% end
% ii
% clear X ii;
%******  Basis Function Generation End******************

mi=0; %minimum iteration
m0=m;   % m0= if Smin becomes minimum in first iteration
I=warp_FFDEPD_3D(m,Moving_Volume,NB,'linear',dxdm,x,y,z,d1,d2,d3);%getting colser to reference
%figure(1), imagesc(I(:,:,2)), colormap(gray)
Smin = inf;  

%sepd=zeros(1,TI);
for iter = 1:Total_Iteration
    [Smin]= update_m_FFDEPD(Smin,I,edge_list_R,iter,d1,d2,d3,NB,dxdm,n,LT,UT);
    %Smin
    I = warp_FFDEPD_3D(m,Moving_Volume,NB,'linear',dxdm,x,y,z,d1,d2,d3);% place I0 at updated position
    %figure(2), imagesc(I(:,:,2)), colormap(gray)
    %iter
end   
if mi==0
    m=m0; %Find m for which epd minimum 
    
else
    m=sm(mi,:);
end
%mi
Registered_Parameters=m;
Registered_Moving_Volume = warp_FFDEPD_3D(m,Moving_Volume,NB,'linear',dxdm,x,y,z,d1,d2,d3);
clear sm m mi Smin m0 I;
%figure(1), imagesc([sum(R0,3)   sum(I,3)]), colormap(gray) 
%V=I;VV=R0;
%Time=toc;
%save FFD54321-DS128-V37VV3 V VV m m0 Smin sm sepd TI mi Time
% figure(2)
% plot(1:TI,sm(:,1),'b',1:TI,sm(:,2),'r',1:TI,sm(:,3),'c',1:TI,sm(:,4),'y',1:TI,sm(:,5),'m',1:TI,sm(:,6),'*',1:TI,sm(:,7),'+',1:TI,sm(:,8),'d',1:TI,sm(:,9),'s',1:TI,sm(:,10),'o',1:TI,sm(:,11),'v',1:TI,sm(:,12),'>')
% figure(3)
% plot(1:TI,sepd(1,:),'g')


