function [Registered_Moving_Volume,Registered_Parameters] = Registration_DPSW_EPD(Moving_V,Total_Iteration,LT,UT,d1,d2,d3,x,y,z,edge_list_R,n,s,dxdmS96) 
% clc;
% close all
% clear;
% warning off all;
global m sm mi
%TI=150; %TI=total iteration
%m=[0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000];
m=[0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000];
sm=zeros(Total_Iteration,3*s^3); %sm= save m, iteration,no. of parameters
%s=2; % s defines no. of parameters, s=2,m=24 s=3, m=81
%LT=0.1;UT=0.9;
%load('Spline_Wavelet.mat') % loading spline wavelet vector SW(1,30001) step size=0.0001, range:0 to 3
%tic


%**************loading reference and moving image volume ********
%load('AffineEPD-DS128-V37VV3Manual.mat')
%  load('Patient 30 original MRI.mat') %RR=R(:,:,11:34); pt 30 ROI %III=II(:,:,10:34); pt56 ROI
%  R0=imresize(CT,.25);
% load('Patient-3CI128128128-MRI.mat')
% R0=V;
% clear V ct_info;

%  load('Patient 56 original MRI.mat')
%  I0=imresize(CT,.25);
% load('Patient-37CI128128128-MRI.mat')
% I0=V;
%%******************** end of data load********************************
%clear V VV;
%clear V ct_info;

% creating edge image of reference volume
%[d1,d2,d3]=size(R0);
%[x,y,z] = meshgrid((0:d2-1),(0:d1-1),(0:d3-1));
%[x,y,z] = meshgrid(1:d2,1:d1,1:d3);

% ER = zeros(d1,d2,d3);
% for i=1:d3
%     ER(:,:,i) = edge(R0(:,:,i),'canny',[LT UT],1.5); % need to select threshold carefully
% end
% [r,c,v] = ind2sub(size(ER),find(ER == 1)); % convert binary edge map into list of pixel positions
% edge_list_R = [c r v];
% [n,~,~] = size(edge_list_R);
% 
% clear r c v ER i;

%********** Basis Function Generation**********

% dxdmS192 = zeros(d1,d2,d3,s^3);
% k = 1;
% for u = 0:s-1
%     for v = 0:s-1
%         for w=0:s-1
%             XW=x*u;
%             XW(find(XW>96)) = XW(find(XW>96))-192*round(XW(find(XW>96))/192); % for 1.5 support and 1.5/2=.75
%             YW=y*v;
%             YW(find(YW>96)) = YW(find(YW>96))-192*round(YW(find(YW>96))/192);
%             ZW=z*w;
%             ZW(find(ZW>96)) = ZW(find(ZW>96))-192*round(ZW(find(ZW>96))/192);
%             dxdmS192(:,:,:,k) = BSpline_wavelet_S192(XW).*BSpline_wavelet_S192(YW).*BSpline_wavelet_S192(ZW);
%             k = k+1;
%         end
%     end
% end

% dxdmS96 = zeros(d1,d2,d3,s^3);
% k = 1;
% for u = 0:s-1
%     for v = 0:s-1
%         for w=0:s-1
%             XW=x*u;
%             XW(find(XW>48)) = XW(find(XW>48))-96*round(XW(find(XW>48))/96); % for 1.5 support and 1.5/2=.75
%             YW=y*v;
%             YW(find(YW>48)) = YW(find(YW>48))-96*round(YW(find(YW>48))/96);
%             ZW=z*w;
%             ZW(find(ZW>48)) = ZW(find(ZW>48))-96*round(ZW(find(ZW>48))/96);
%             dxdmS96(:,:,:,k) = BSpline_wavelet_S96(XW).*BSpline_wavelet_S96(YW).*BSpline_wavelet_S96(ZW);
%             k = k+1;
%             clear  XW YW ZW;
%         end
%     end
% end
% 
% % dxdm = zeros(d1,d2,d3,s^3);
% % for k=1:s^3 
% % dxdm(:,:,:,k)= squeeze(dxdmS192(:,:,:,k))+squeeze(dxdmS96(:,:,:,k));  
% % end
% 
% clear k u v w;


%******End******************

mi=0; %minimum iteration
m0=m;   % m0= if Smin becomes minimum in first iteration
I=warp_DPSWEPD_3D(m,Moving_V,s,'linear',dxdmS96,x,y,z,d1,d2,d3);%getting colser to reference
Smin = inf;  

%sepd=zeros(1,TI);
for iter = 1:Total_Iteration
    [Smin]= update_m_DPSWEPD(Smin,I,edge_list_R,iter,d1,d2,d3,s,dxdmS96,LT,UT,n);
    %Smin
    I = warp_DPSWEPD_3D(m,Moving_V,s,'linear',dxdmS96,x,y,z,d1,d2,d3);% place I0 at updated position
    %iter
end   
if mi==0
    m=m0; %Find m for which epd minimum 
else
    m=sm(mi,:);
end
%mi
Registered_Parameters=m;
Registered_Moving_Volume = warp_DPSWEPD_3D(m,Moving_V,s,'linear',dxdmS96,x,y,z,d1,d2,d3);
clear sm m mi Smin m0 I;
% figure(1), imagesc([sum(R0,3)   sum(I,3)]), colormap(gray) 
% V=I;VV=R0;
%Time=toc;
%save DPSWEPD-DS128-V37VV3 V VV m m0 Smin sm sepd TI mi Time 
% figure(2)
% plot(1:TI,sm(:,1),'b',1:TI,sm(:,2),'r',1:TI,sm(:,3),'c',1:TI,sm(:,4),'y',1:TI,sm(:,5),'m',1:TI,sm(:,6),'*',1:TI,sm(:,7),'+',1:TI,sm(:,8),'d',1:TI,sm(:,9),'s',1:TI,sm(:,10),'o',1:TI,sm(:,11),'v',1:TI,sm(:,12),'>')
% figure(3)
% plot(1:TI,sepd(1,:),'g')


