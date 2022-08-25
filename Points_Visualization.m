clc
clear all

% Fixed data
load('Patient 3CI muscle data.mat')
[Xf,Yf,Zf] = get_ptCloud(muscle_data);

% Moving data
load('Patient 13CI muscle data.mat')
[Xm,Ym,Zm] = get_ptCloud(muscle_data);
figure(1); plot_ptCloud(Xf,Yf,Zf,Xm,Ym,Zm)

% loading registered parameters
load('Registration_Parameters_3F_13M_xCornerO.mat')


% initialization
d1=128;d2=128;d3=128;
[x,y,z] = meshgrid(1:d2,1:d1,1:d3);
dxdmS96 = calculate_dxdm_DPSW(2,d1,d2,d3,x,y,z);

j5=5;
NB5=((2^-j5)*d1+1)*1;
dxdm5 = calculate_dxdm_FFD(d1,d2,d3,NB5,j5);
j4=4; NB4=((2^-j4)*d1+1)*1;
dxdm4 = calculate_dxdm_FFD(d1,d2,d3,NB4,j4);
j3=3; NB3=((2^-j3)*d1+1)*1;
dxdm3 = calculate_dxdm_FFD(d1,d2,d3,NB3,j3);
j2=2; NB2=((2^-j2)*d1+1)*1;
dxdm2 = calculate_dxdm_FFD(d1,d2,d3,NB2,j2);
j1=1; NB1=((2^-j1)*d1+1)*1;
dxdm1 = calculate_dxdm_FFD(d1,d2,d3,NB1,j1);

% Transfereing points
[Xf_TH,Yf_TH,Zf_TH] = PointTrans_Affine_EPD(R_Parameters_AffineEPD,Xf,Yf,Zf);
figure(2); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)

[Xf_TH,Yf_TH,Zf_TH] = PointTrans_DPSWEPD(R_Parameters_DPSWEPD,Xf_TH,Yf_TH,Zf_TH,dxdmS96);
figure(3); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)

[Xf_TH,Yf_TH,Zf_TH] = PointTrans_FFD_EPD(R_Parameters_FFD5_EPD,Xf_TH,Yf_TH,Zf_TH,dxdm5,NB5);
figure(4); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)

[Xf_TH,Yf_TH,Zf_TH] = PointTrans_FFD_EPD(R_Parameters_FFD4_EPD,Xf_TH,Yf_TH,Zf_TH,dxdm4,NB4);
figure(5); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)

[Xf_TH,Yf_TH,Zf_TH] = PointTrans_FFD_EPD(R_Parameters_FFD3_EPD,Xf_TH,Yf_TH,Zf_TH,dxdm3,NB3);
figure(6); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)

[Xf_TH,Yf_TH,Zf_TH] = PointTrans_FFD_EPD(R_Parameters_FFD2_EPD,Xf_TH,Yf_TH,Zf_TH,dxdm2,NB2);
figure(7); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)

[Xf_TH,Yf_TH,Zf_TH] = PointTrans_FFD_EPD(R_Parameters_FFD1_EPD,Xf_TH,Yf_TH,Zf_TH,dxdm1,NB1);
figure(8); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)

[Xf_TH,Yf_TH,Zf_TH] = PointTrans_Demons(Displacements_Hybrid,Xf_TH,Yf_TH,Zf_TH);
figure(9); plot_ptCloud(Xf_TH,Yf_TH,Zf_TH,Xm,Ym,Zm)

[Xf_TD,Yf_TD,Zf_TD] = PointTrans_Demons(Displacements_Demons,Xf,Yf,Zf);
figure(10); plot_ptCloud(Xf_TD,Yf_TD,Zf_TD,Xm,Ym,Zm)

save transfered_points_3F_13M Xf_TD Yf_TD Zf_TD Xf_TH Yf_TH Zf_TH
