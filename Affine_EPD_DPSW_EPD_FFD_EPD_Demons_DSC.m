%Pleas contact SUMAN(mobile:0424364536) before use this computer
% Building-15, Room-153

% This code calculate registration framework and DSC for single fixed image volume in respect of multiple Moving image volumes.  
clc;
close all
clear;
warning off all;
tic

% Loading Fixed & Moving volumes and muscles data

load('Patient-61CI128128128-MRI.mat')
load('Patient 61CI muscle data.mat')
Fixed_Volume=V; % Fixed volume or Reference volume
F_muscle_data=muscle_data; %[Xf,Yf,Zf] = get_ptCloud(F_muscle_data);
clear V ct_info muscle_data;


load('Patient-3CI128128128-MRI.mat')
load('Patient 3CI muscle data.mat')
Moving_Volume{1}=V;M_muscle_data{1}=muscle_data; clear V ct_info muscle_data; % moving Volume
load('Patient-13CI128128128-MRI.mat')
load('Patient 13CI muscle data.mat')
Moving_Volume{2}=V;M_muscle_data{2}=muscle_data;clear V ct_info muscle_data;
load('Patient-16CI128128128-MRI.mat')
load('Patient 16CI muscle data.mat')
Moving_Volume{3}=V;M_muscle_data{3}=muscle_data;clear V ct_info muscle_data;
load('Patient-17CI128128128-MRI.mat')
load('Patient 17CI muscle data.mat')
Moving_Volume{4}=V;M_muscle_data{4}=muscle_data;clear V ct_info muscle_data;
load('Patient-30CI128128128-MRI.mat')
load('Patient 30CI muscle data.mat')
Moving_Volume{5}=V;M_muscle_data{5}=muscle_data;clear V ct_info muscle_data;
load('Patient-33CI128128128-MRI.mat')
load('Patient 33CI muscle data.mat')
Moving_Volume{6}=V;M_muscle_data{6}=muscle_data;clear V ct_info muscle_data; 
load('Patient-34CI128128128-MRI.mat')
load('Patient 34CI muscle data.mat')
Moving_Volume{7}=V;M_muscle_data{7}=muscle_data;clear V ct_info muscle_data;
load('Patient-37CI128128128-MRI.mat')
load('Patient 37CI muscle data.mat')
Moving_Volume{8}=V;M_muscle_data{8}=muscle_data;clear V ct_info muscle_data;
load('Patient-45CI128128128-MRI.mat')
load('Patient 45CI muscle data.mat')
Moving_Volume{9}=V;M_muscle_data{9}=muscle_data;clear V ct_info muscle_data;
load('Patient-56CI128128128-MRI.mat')
load('Patient 56CI muscle data.mat')
Moving_Volume{10}=V;M_muscle_data{10}=muscle_data;clear V ct_info muscle_data; 

% % % Creating matrix for holding DSC values
TH_DT6502 = zeros(10,6); % TH_DT=tabke hybrid dual threshold, 10 moving, 6 muscles, DSC tables
TH_DT102 = zeros(10,6);
TH_ST01=zeros(10,6);
TH_ST2=zeros(10,6);
TH_ST3=zeros(10,6);
TH_ST5 = zeros(10,6); % ST= single threshold
TH_ST65 = zeros(10,6);
TH_ST75=zeros(10,6);
TH_ST1 = zeros(10,6);
TH_ST5A1CA= zeros(10,6);
TH_ST5A1C= zeros(10,6);
TD_T02 = zeros(10,6);
TD_T005 = zeros(10,6);
TD_T5 = zeros(10,6);

%Initialiazation for whole framework
[d1,d2,d3]=size(Fixed_Volume);
[x,y,z] = meshgrid(1:d2,1:d1,1:d3); %  coordinates as origin top left corner of the volume; intrinsic or spatial coordinate system
[xOC,yOC,zOC] = meshgrid((0:d2-1)-d2/2,(0:d1-1)-d1/2,(0:d3-1)-d3/2); % coordinates as origin center of the volume; intrinsic or spatial coordinate system
TI150=150; %TI=total iteration
TI100=100; TI50=50;
LT09=0.09;UT7=0.7;% Lower thresold
[Edge_list097,Number_EdgePoints097] = Calculate_EdgeList(Fixed_Volume,d1,d2,d3,LT09,UT7);
LT08=0.08;%UT7=0.7;
[Edge_list087,Number_EdgePoints087] = Calculate_EdgeList(Fixed_Volume,d1,d2,d3,LT08,UT7);
LT04=0.04;UT4=0.4;
[Edge_list044,Number_EdgePoints044] = Calculate_EdgeList(Fixed_Volume,d1,d2,d3,LT04,UT4);
LT01=0.01;UT2=0.2;
[Edge_list012,Number_EdgePoints012] = Calculate_EdgeList(Fixed_Volume,d1,d2,d3,LT01,UT2);
UT1=0.1;
[Edge_list011,Number_EdgePoints011] = Calculate_EdgeList(Fixed_Volume,d1,d2,d3,LT01,UT1);



% Initialization for DPSW5(5th resolution) Registration
BFD=2; % Basis function determinant, 2 means 8 basis functions & 8*3=24 parameters; 3 means 27 basis functions & 27*3=81 parameters
dxdmS96 = calculate_dxdm_DPSW_S96(BFD,d1,d2,d3,x,y,z); % 5th resolution means support=96 of wavelet, 3 original support, 3*2^5=96

% Initialization for DPSW4(4th resolution) Registration
dxdmS48 = calculate_dxdm_DPSW_S48(BFD,d1,d2,d3,x,y,z);
    
% Initialization for FFD5 Registration 
j5=5;% original spline support=4, jj=6 means support=256, 256/4=64=2^6, Resolution Level 
NB5=((2^-j5)*d1+1)*1; % number of basis functions for 5th resolution FFD
dxdm5 = calculate_dxdm_FFD(d1,d2,d3,NB5,j5);

% Initialization for FFD4 Registration 
j4=4; NB4=((2^-j4)*d1+1)*1;
dxdm4 = calculate_dxdm_FFD(d1,d2,d3,NB4,j4);

% Initialization for FFD3 Registration 
j3=3; NB3=((2^-j3)*d1+1)*1;
dxdm3 = calculate_dxdm_FFD(d1,d2,d3,NB3,j3);

% Initialization for FFD2 Registration 
j2=2; NB2=((2^-j2)*d1+1)*1;
dxdm2 = calculate_dxdm_FFD(d1,d2,d3,NB2,j2);

% Initialization for FFD1 EPD Registration 
j1=1; NB1=((2^-j1)*d1+1)*1;
dxdm1 = calculate_dxdm_FFD(d1,d2,d3,NB1,j1);


for MVI=1:10 % MVI=moving volume index
    %tic
    %[Xm,Ym,Zm] = get_ptCloud(M_muscle_data{moving_volume_index});
    
%     % Corner origin to Center Origin Conversion
%     [Xf,Yf,Zf]=ptCloud_OriTra(Xf,Yf,Zf,d1,d2,d3);
%     [Xm,Ym,Zm]=ptCloud_OriTra(Xm,Ym,Zm,d1,d2,d3);  
%     figure(1); plot_ptCloud(Xf,Yf,Zf,Xm,Ym,Zm)
    
    % Window Section
    %[RV_WindowMI,RP_WindowMI] = Registration_Window_MI(Fixed_Volume,,xOC,yOC,zOC);
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_ST(R_P_WindowMI,Xf,Yf,Zf);
%     figure(2); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    %[RV_WindowEPD,RP_WindowEPD] = Registration_Window_EPD(Fixed_Volume,RV_WindowMI,d1,d2,d3,LT1,UT9,xOC,yOC,zOC);
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_ST(R_P_WindowEPD,Xf_TH,Yf_TH,Zf_TH);
%     figure(3); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    %clear  RV_WindowMI;
    
    %1st Section
    [RV_Affine19,RP_Affine19] = Registration_Affine_EPD(Moving_Volume{MVI},TI150,LT09,UT7,d1,d2,d3,xOC,yOC,zOC,Edge_list097,Number_EdgePoints097);
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_Affine_EPD(R_P_Affine19,Xf_TH,Yf_TH,Zf_TH);
%     figure(4); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    %clear  RV_WindowEPD;
    
    [DisH_A1CA,RV_HDA1CA] = imregdemons(RV_Affine19,Fixed_Volume);
    
    [RV_DPSW5_19,RP_DPSW5_19] = Registration_DPSW_EPD(RV_Affine19,TI150,LT09,UT7,d1,d2,d3,xOC,yOC,zOC,Edge_list097,Number_EdgePoints097,BFD,dxdmS96);
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_DPSWEPD(R_P_DPSW5_19,Xf_TH,Yf_TH,Zf_TH,dxdmS96,d1,d2,d3);
%     figure(5); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    clear RV_Affine19;
%     [R_M_V_DPSW4_19,R_P_DPSW4_19] = Registration_DPSW_EPD(R_M_V_DPSW5_19,TI150,LT1,UT9,d1,d2,d3,xOC,yOC,zOC,Edge_list19,Number_EdgePoints19,BFD,dxdmS48);
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_DPSWEPD(R_P_DPSW4_19,Xf_TH,Yf_TH,Zf_TH,dxdmS48,d1,d2,d3);
%     figure(6); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    %clear R_M_V_DPSW5_19;   
    [RV_FFD5,RP_FFD5] = Registration_FFD_EPD(RV_DPSW5_19,TI100,LT09,UT7,d1,d2,d3,xOC,yOC,zOC,Edge_list097,Number_EdgePoints097,NB5,dxdm5);
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_FFD_EPD(R_P_FFD5,Xf_TH,Yf_TH,Zf_TH,dxdm5,NB5,d1,d2,d3);
%     figure(6); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    clear RV_DPSW5_19;
    
    [DisH_A1C,RV_HDA1C] = imregdemons(RV_FFD5,Fixed_Volume);
    
    %2nd Section
    %[RV_Affine087,RP_Affine087] = Registration_Affine_EPD(RV_FFD5,TI50,LT08,UT7,d1,d2,d3,xOC,yOC,zOC,Edge_list087,Number_EdgePoints087);
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_Affine_EPD(R_P_Affine087,Xf_TH,Yf_TH,Zf_TH);
%     figure(7); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    %clear  RV_FFD5;
    %[RV_DPSW5_087,RP_DPSW5_087] = Registration_DPSW_EPD(RV_Affine087,TI50,LT08,UT7,d1,d2,d3,xOC,yOC,zOC,Edge_list087,Number_EdgePoints087,BFD,dxdmS96); 
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_DPSWEPD(R_P_DPSW5_087,Xf_TH,Yf_TH,Zf_TH,dxdmS96,d1,d2,d3);
%     figure(8); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    %clear RV_Affine087;
    %[RV_DPSW4_087,RP_DPSW4_087] = Registration_DPSW_EPD(RV_DPSW5_087,TI50,LT08,UT7,d1,d2,d3,xOC,yOC,zOC,Edge_list087,Number_EdgePoints087,BFD,dxdmS48); 
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_DPSWEPD(R_P_DPSW4_087,Xf_TH,Yf_TH,Zf_TH,dxdmS48,d1,d2,d3);
%     figure(9); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    %clear RV_DPSW5_087;   
    [RV_FFD4,RP_FFD4] = Registration_FFD_EPD(RV_FFD5,TI100,LT08,UT7,d1,d2,d3,xOC,yOC,zOC,Edge_list087,Number_EdgePoints087,NB4,dxdm4);
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_FFD_EPD(R_P_FFD4,Xf_TH,Yf_TH,Zf_TH,dxdm4,NB4,d1,d2,d3);
%     figure(10); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    clear RV_FFD5;
    
    %[DisH_A2S,~] = imregdemons(RV_FFD4,Fixed_Volume);
    
    %3rd Section
    %[RV_Affine044,RP_Affine044] = Registration_Affine_EPD(RV_FFD4,TI50,LT04,UT4,d1,d2,d3,xOC,yOC,zOC,Edge_list044,Number_EdgePoints044);
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_Affine_EPD(R_P_Affine044,Xf_TH,Yf_TH,Zf_TH);
%     figure(11); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    %clear  RV_FFD4;
    %[RV_DPSW5_044,RP_DPSW5_044] = Registration_DPSW_EPD(RV_Affine044,TI50,LT04,UT4,d1,d2,d3,xOC,yOC,zOC,Edge_list044,Number_EdgePoints044,BFD,dxdmS96); 
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_DPSWEPD(R_P_DPSW5_044,Xf_TH,Yf_TH,Zf_TH,dxdmS96,d1,d2,d3);
%     figure(12); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    %clear RV_Affine044;
    %[RV_DPSW4_044,RP_DPSW4_044] = Registration_DPSW_EPD(RV_DPSW5_044,TI50,LT04,UT4,d1,d2,d3,xOC,yOC,zOC,Edge_list044,Number_EdgePoints044,BFD,dxdmS48); 
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_DPSWEPD(R_P_DPSW4_044,Xf_TH,Yf_TH,Zf_TH,dxdmS48,d1,d2,d3);
%     figure(13); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    %clear RV_DPSW5_044;   
    [RV_FFD3,RP_FFD3] = Registration_FFD_EPD(RV_FFD4,TI100,LT04,UT4,d1,d2,d3,xOC,yOC,zOC,Edge_list044,Number_EdgePoints044,NB3,dxdm3);
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_FFD_EPD(R_P_FFD3,Xf_TH,Yf_TH,Zf_TH,dxdm3,NB3,d1,d2,d3);
%     figure(14); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    clear RV_FFD4;
    
    %[DisH_A3S,~] = imregdemons(RV_FFD3,Fixed_Volume);
    
    %4th Section
    %[R_M_V_Affine012,R_P_Affine012] = Registration_Affine_EPD(R_M_V_FFD3,TI50,LT01,UT2,d1,d2,d3,xOC,yOC,zOC,Edge_list012,Number_EdgePoints012);
    %clear  R_M_V_FFD3;
    %[R_M_V_DPSW5_012,R_P_DPSW5_012] = Registration_DPSW_EPD(R_M_V_Affine012,TI50,LT01,UT2,d1,d2,d3,xOC,yOC,zOC,Edge_list012,Number_EdgePoints012,BFD,dxdmS96); 
    %clear R_M_V_Affine012;
    %[R_M_V_DPSW4_012,R_P_DPSW4_012] = Registration_DPSW_EPD(R_M_V_DPSW5_012,TI50,LT01,UT2,d1,d2,d3,xOC,yOC,zOC,Edge_list012,Number_EdgePoints012,BFD,dxdmS48); 
    %clear R_M_V_DPSW5_012;   
    [RV_FFD2,RP_FFD2] = Registration_FFD_EPD(RV_FFD3,TI100,LT01,UT2,d1,d2,d3,xOC,yOC,zOC,Edge_list012,Number_EdgePoints012,NB2,dxdm2);
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_FFD_EPD(R_P_FFD2,Xf_TH,Yf_TH,Zf_TH,dxdm2,NB2,d1,d2,d3);
%     figure(15); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    clear RV_FFD3;
    
    %5th Section
    %[R_M_V_Affine011,R_P_Affine011] = Registration_Affine_EPD(R_M_V_FFD2,TI50,LT01,UT1,d1,d2,d3,xOC,yOC,zOC,Edge_list011,Number_EdgePoints011);
    %clear  R_M_V_FFD2;
    %[R_M_V_DPSW5_011,R_P_DPSW5_011] = Registration_DPSW_EPD(R_M_V_Affine011,TI50,LT01,UT1,d1,d2,d3,xOC,yOC,zOC,Edge_list011,Number_EdgePoints011,BFD,dxdmS96); 
    %clear R_M_V_Affine011;
    %[R_M_V_DPSW4_011,R_P_DPSW4_011] = Registration_DPSW_EPD(R_M_V_DPSW5_011,TI50,LT01,UT1,d1,d2,d3,xOC,yOC,zOC,Edge_list011,Number_EdgePoints011,BFD,dxdmS48); 
    %clear R_M_V_DPSW5_011;   
    [RV_FFD1,RP_FFD1] = Registration_FFD_EPD(RV_FFD2,TI100,LT01,UT1,d1,d2,d3,xOC,yOC,zOC,Edge_list011,Number_EdgePoints011,NB1,dxdm1);
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_FFD_EPD(R_P_FFD1,Xf_TH,Yf_TH,Zf_TH,dxdm1,NB1,d1,d2,d3);
%     figure(16); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    clear RV_FFD2;
    
    [DisH,RV_HDAF] = imregdemons(RV_FFD1,Fixed_Volume);
%     [Xf_TH,Yf_TH,Zf_TH] = PointTrans_Demons(Displacements_Hybrid,Xf_TH,Yf_TH,Zf_TH,d1,d2,d3);
%     figure(17); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    %clear RV_FFD1;
    [DisD,RV_D] = imregdemons(Moving_Volume{MVI},Fixed_Volume);
%     [Xf_TD,Yf_TD,Zf_TD] = PointTrans_Demons(Displacements_Demons,Xf,Yf,Zf,d1,d2,d3);
%     figure(18); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)
    
%      % Center origin to Corner Origin Conversion
%     [Xm,Ym,Zm]=ptCloud_OriTra_Cen2Con(Xm,Ym,Zm,d1,d2,d3);
%     [Xf_TH,Yf_TH,Zf_TH]=ptCloud_OriTra_Cen2Con(Xf_TH,Yf_TH,Zf_TH,d1,d2,d3);
%     [Xf_TD,Yf_TD,Zf_TD]=ptCloud_OriTra_Cen2Con(Xf_TD,Yf_TD,Zf_TD,d1,d2,d3);
    
%     [DSC_Table_Hybrid(moving_volume_index,:),DSC_Table_Demons(moving_volume_index,:)] = DSC_Cal_Points_Warp(Xm,Ym,Zm,Xf_TH,Yf_TH,Zf_TH,Xf_TD,Yf_TD,Zf_TD,d1,d2,d3); 
    
     %save 3F_16M RV_Affine19 RV_DPSW5_19 RV_FFD5 RV_FFD4 RV_FFD3 RV_FFD2 RV_FFD1 RV_D RV_HDAF RV_HDA1C RV_HDA1CA
%     toc

%     [TH_DT6502(MVI,:),TH_DT102(MVI,:),TD_T02(MVI,:)] = Vol_LabelT_Warp_DSC_Cal_DT(F_muscle_data,M_muscle_data{MVI},d1,d2,d3,xOC,yOC,zOC,RP_WindowMI,RP_WindowEPD,RP_Affine19,RP_FFD5,RP_Affine087,RP_FFD4,RP_Affine044,RP_FFD3,RP_FFD2,RP_FFD1,DisH,DisD,dxdm5,dxdm4,dxdm3,dxdm2,dxdm1,NB5,NB4,NB3,NB2,NB1);   
%     [TH_ST01(MVI,:),TH_ST2(MVI,:),TH_ST3(MVI,:),TH_ST5(MVI,:),TH_ST65(MVI,:),TH_ST75(MVI,:),TH_ST1(MVI,:),TH_ST5A2S(MVI,:),TD_T005(MVI,:),TD_T5(MVI,:)] = Vol_LabelT_Warp_DSC_Cal_ST(F_muscle_data,M_muscle_data{MVI},d1,d2,d3,xOC,yOC,zOC,RP_WindowMI,RP_WindowEPD,RP_Affine19,RP_FFD5,RP_Affine087,RP_FFD4,RP_Affine044,RP_FFD3,RP_FFD2,RP_FFD1,DisH,DisD,DisH_A2S,dxdm5,dxdm4,dxdm3,dxdm2,dxdm1,NB5,NB4,NB3,NB2,NB1); 
%     clear RP_WindowMI RP_WindowEPD RP_Affine19 RP_FFD5 RP_Affine087 RP_FFD4 RP_Affine044 RP_FFD3 RP_FFD2 RP_FFD1 DisD DisH DisH_A2S;
    
    [TH_DT6502(MVI,:),TH_DT102(MVI,:),TD_T02(MVI,:)] = Vol_LabelT_Warp_DSC_Cal_DT(F_muscle_data,M_muscle_data{MVI},d1,d2,d3,xOC,yOC,zOC,RP_Affine19,RP_DPSW5_19,RP_FFD5,RP_FFD4,RP_FFD3,RP_FFD2,RP_FFD1,DisH,DisD,dxdmS96,dxdm5,dxdm4,dxdm3,dxdm2,dxdm1,BFD,NB5,NB4,NB3,NB2,NB1);   
    [TH_ST01(MVI,:),TH_ST2(MVI,:),TH_ST3(MVI,:),TH_ST5(MVI,:),TH_ST65(MVI,:),TH_ST75(MVI,:),TH_ST1(MVI,:),TH_ST5A1CA(MVI,:),TH_ST5A1C(MVI,:),TD_T005(MVI,:),TD_T5(MVI,:)] = Vol_LabelT_Warp_DSC_Cal_ST(F_muscle_data,M_muscle_data{MVI},d1,d2,d3,xOC,yOC,zOC,RP_Affine19,RP_DPSW5_19,RP_FFD5,RP_FFD4,RP_FFD3,RP_FFD2,RP_FFD1,DisH,DisD,DisH_A1CA,DisH_A1C,dxdmS96,dxdm5,dxdm4,dxdm3,dxdm2,dxdm1,BFD,NB5,NB4,NB3,NB2,NB1); 
    clear RP_Affine19 RP_DPSW5_19 RP_FFD5 RP_FFD4 RP_FFD3 RP_FFD2 RP_FFD1 DisD DisH DisH_A1CA DisH_A1C;
    %save P61F56M097CS_SFS TH_DT6502 TH_DT102 TH_ST01 TH_ST2 TH_ST3 TH_ST5 TH_ST65 TH_ST75 TH_ST1 TH_ST5A1CA TH_ST5A1C TD_T005 TD_T5 TD_T02
    MVI
end
Time=toc;
save P61F TH_DT6502 TH_DT102 TH_ST01 TH_ST2 TH_ST3 TH_ST5 TH_ST65 TH_ST75 TH_ST1 TH_ST5A1CA TH_ST5A1C TD_T005 TD_T5 TD_T02 Time  













