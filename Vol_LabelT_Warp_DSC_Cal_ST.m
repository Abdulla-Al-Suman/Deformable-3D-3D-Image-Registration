%function [VH01,VH2,VH3,VH5,VH65,VH75,VH1,VH5A2S,VD005,VD5] = Vol_LabelT_Warp_DSC_Cal_ST(F_muscle_data,M_muscle_data,d1,d2,d3,xOC,yOC,zOC,RP_WindowMI,RP_WindowEPD,RP_Affine19,RP_FFD5,RP_Affine087,RP_FFD4,RP_Affine044,RP_FFD3,RP_FFD2,RP_FFD1,DisH,DisD,DisH_A2S,dxdm5,dxdm4,dxdm3,dxdm2,dxdm1,NB5,NB4,NB3,NB2,NB1)
function [VH01,VH2,VH3,VH5,VH65,VH75,VH1,VH5A1CA,VH5A1C,VD005,VD5] = Vol_LabelT_Warp_DSC_Cal_ST(F_muscle_data,M_muscle_data,d1,d2,d3,xOC,yOC,zOC,RP_Affine19,RP_DPSW5_19,RP_FFD5,RP_FFD4,RP_FFD3,RP_FFD2,RP_FFD1,DisH,DisD,DisH_A1CA,DisH_A1C,dxdmS96,dxdm5,dxdm4,dxdm3,dxdm2,dxdm1,BFD,NB5,NB4,NB3,NB2,NB1)
%VH5=vector hybrid
% transfering Labels using volume template process

% Creating vector for holding DSC values
DSC_VH01 = zeros(1,6);
DSC_VH2 = zeros(1,6);
DSC_VH3 = zeros(1,6);
DSC_VH5 = zeros(1,6);
DSC_VH65 = zeros(1,6);
DSC_VH75 = zeros(1,6);
DSC_VH1 = zeros(1,6);
DSC_VH5A1CA = zeros(1,6);
DSC_VH5A1C = zeros(1,6);
DSC_VD005 = zeros(1,6);
DSC_VD5 = zeros(1,6);


muscle = [3 16 22]; % muscles numbers
for MN = 1:length(muscle) % MN=muscle number
    for SN = 1:2 %SNA=1(Left)= side number,
        
        % Fixed image mask generation on ROI
        Fixed_Mask=zeros(d1,d2,d3);
        for FN = 40:70  
            X = F_muscle_data{muscle(MN)}{SN}{FN}.x; 
            Y = F_muscle_data{muscle(MN)}{SN}{FN}.y;
            Fixed_Mask(:,:,FN) = poly2mask(X,Y,d2,d1);
            clear FN X Y;
        end
        
        % Moving image banary volume generation
        Moving=zeros(d1,d2,d3);
        for FN = 1:d3  % 1:O for full data, FN= frame number, for making banary volume
            X = M_muscle_data{muscle(MN)}{SN}{FN}.x;  % muscle selection, annotation software save x as 1 to N
            Y = M_muscle_data{muscle(MN)}{SN}{FN}.y;  % 1= left, 2=right in the original annotation software
            Moving(:,:,FN) = poly2mask(X,Y,d2,d1);
            clear FN X Y;
        end  
        
        % data transfer for  only Demon
        IODemons = imwarp(Moving,DisD);
        IIODemons005 = ImageThresholding(IODemons,0.005,d1,d2,d3);
        IIODemons5 = ImageThresholding(IODemons,0.5,d1,d2,d3);
        clear IODemons;
        
        % START Window Search
        % data transfer for Window MI transformation
%         IW_MI = warp_Window_ST_3D(Moving,RP_WindowMI,xOC,yOC,zOC); 
%         clear Moving;        
%        % data transfer for Window EPD transformation
%         IW_EPD = warp_Window_ST_3D(IW_MI,RP_WindowEPD,xOC,yOC,zOC);
%         clear IW_MI;
        
         % START 1st Section
        % data transfer for affine transformation
        IAffine19 = warp_AffineEPD_3D(Moving,RP_Affine19,xOC,yOC,zOC);
        clear Moving;   
        
        IHybridA1CA = imwarp(IAffine19,DisH_A1CA);
        IHybrid5A1CA = ImageThresholding(IHybridA1CA,.5,d1,d2,d3);
        clear IHybridA1CA;
        
        % data transfer for DPSW5 transformation
        IDPSW5_19 = warp_DPSWEPD_3D(RP_DPSW5_19,IAffine19,BFD,'linear',dxdmS96,xOC,yOC,zOC,d1,d2,d3);
        clear IAffine19;
%         % data transfer for DPSW4 transformation
%         IDPSW4_19 = warp_DPSWEPD_3D(R_P_DPSW4_19,IDPSW5_19,BFD,'linear',dxdmS48,xOC,yOC,zOC,d1,d2,d3);
%         clear IDPSW5_19;
        % data transfer for FFD5 transformation
        IFFD5=warp_FFDEPD_3D(RP_FFD5,IDPSW5_19,NB5,'linear',dxdm5,xOC,yOC,zOC,d1,d2,d3);
        clear IDPSW5_19;
        
        IHybridA1C = imwarp(IFFD5,DisH_A1C);
        IHybrid5A1C = ImageThresholding(IHybridA1C,.5,d1,d2,d3);
        clear IHybridA1C;
        
        % START 2nd Section
       % data transfer for affine transformation
%         IAffine087 = warp_AffineEPD_3D(IFFD5,RP_Affine087,xOC,yOC,zOC);
%         clear IFFD5;      
%         % data transfer for DPSW5 transformation
%         IDPSW5_087 = warp_DPSWEPD_3D(RP_DPSW5_087,IAffine087,BFD,'linear',dxdmS96,xOC,yOC,zOC,d1,d2,d3);
%         clear IAffine087;
        % data transfer for DPSW4 transformation
%         IDPSW4_087 = warp_DPSWEPD_3D(RP_DPSW4_087,IDPSW5_087,BFD,'linear',dxdmS48,xOC,yOC,zOC,d1,d2,d3);
%         clear IDPSW5_087;
        % data transfer for FFD54 transformation
        IFFD4=warp_FFDEPD_3D(RP_FFD4,IFFD5,NB4,'linear',dxdm4,xOC,yOC,zOC,d1,d2,d3);
        clear IFFD5;
        
        
%         IHybridA2S = imwarp(IFFD4,DisH_A2S);
%         IHybrid5A2S = ImageThresholding(IHybridA2S,.5,d1,d2,d3);
%         clear IHybridA2S;
        
        % START 3rd Section
        % data transfer for affine transformation
%         IAffine044 = warp_AffineEPD_3D(IFFD4,RP_Affine044,xOC,yOC,zOC);
%         clear IFFD4;      
%         % data transfer for DPSW5 transformation
%         IDPSW5_044 = warp_DPSWEPD_3D(RP_DPSW5_044,IAffine044,BFD,'linear',dxdmS96,xOC,yOC,zOC,d1,d2,d3);
%         clear IAffine044;
        % data transfer for DPSW4 transformation
%         IDPSW4_044 = warp_DPSWEPD_3D(RP_DPSW4_044,IDPSW5_044,BFD,'linear',dxdmS48,xOC,yOC,zOC,d1,d2,d3);
%         clear IDPSW5_044;
        % data transfer for FFD543 transformation
        IFFD3=warp_FFDEPD_3D(RP_FFD3,IFFD4,NB3,'linear',dxdm3,xOC,yOC,zOC,d1,d2,d3);
        clear IFFD4;
        
%         IHybridA3S = imwarp(IFFD3,DisH_A3S);
%         IHybrid5A3S = ImageThresholding(IHybridA3S,.5,d1,d2,d3);
%         clear IHybridA3S;
        
         % START 4th Section
        % data transfer for affine transformation
%         IAffine012 = warp_AffineEPD_3D(IFFD3,R_P_Affine012,xOC,yOC,zOC);
%         clear IFFD3;      
%         % data transfer for DPSW5 transformation
%         IDPSW5_012 = warp_DPSWEPD_3D(R_P_DPSW5_012,IAffine012,BFD,'linear',dxdmS96,xOC,yOC,zOC,d1,d2,d3);
%         clear IAffine012;
%         % data transfer for DPSW4 transformation
%         IDPSW4_012 = warp_DPSWEPD_3D(R_P_DPSW4_012,IDPSW5_012,BFD,'linear',dxdmS48,xOC,yOC,zOC,d1,d2,d3);
%         clear IDPSW5_012;
        % data transfer for FFD5432 transformation
        IFFD2=warp_FFDEPD_3D(RP_FFD2,IFFD3,NB2,'linear',dxdm2,xOC,yOC,zOC,d1,d2,d3);
        clear IFFD3;
        
        % START 5th Section
        % data transfer for affine transformation
%         IAffine011 = warp_AffineEPD_3D(IFFD2,R_P_Affine011,xOC,yOC,zOC);
%         clear IFFD2;      
%         % data transfer for DPSW5 transformation
%         IDPSW5_011 = warp_DPSWEPD_3D(R_P_DPSW5_011,IAffine011,BFD,'linear',dxdmS96,xOC,yOC,zOC,d1,d2,d3);
%         clear IAffine011;
%         % data transfer for DPSW4 transformation
%         IDPSW4_011 = warp_DPSWEPD_3D(R_P_DPSW4_011,IDPSW5_011,BFD,'linear',dxdmS48,xOC,yOC,zOC,d1,d2,d3);
%         clear IDPSW5_011;
        % data transfer for FFD54321 transformation
        IFFD1=warp_FFDEPD_3D(RP_FFD1,IFFD2,NB1,'linear',dxdm1,xOC,yOC,zOC,d1,d2,d3);
        clear IFFD2;
        
        % data transfer for Hybrid Demon
        IHybrid = imwarp(IFFD1,DisH);
        clear IFFD1;
        IHybrid01 = ImageThresholding(IHybrid,.1,d1,d2,d3);
        IHybrid2 = ImageThresholding(IHybrid,.2,d1,d2,d3);
        IHybrid3 = ImageThresholding(IHybrid,.3,d1,d2,d3);
        IHybrid5 = ImageThresholding(IHybrid,.5,d1,d2,d3);
        IHybrid65 = ImageThresholding(IHybrid,.65,d1,d2,d3);
        IHybrid75 = ImageThresholding(IHybrid,.75,d1,d2,d3);
        IHybrid1 = ImageThresholding(IHybrid,1,d1,d2,d3);
        clear IHybrid;
   
        
        % making ROI warped moving Image
        Warped_Moving_Mask_Hybrid01=zeros(d1,d2,d3);
        Warped_Moving_Mask_Hybrid2=zeros(d1,d2,d3);
        Warped_Moving_Mask_Hybrid3=zeros(d1,d2,d3);
        Warped_Moving_Mask_Hybrid5=zeros(d1,d2,d3);
        Warped_Moving_Mask_Hybrid65=zeros(d1,d2,d3);
        Warped_Moving_Mask_Hybrid75=zeros(d1,d2,d3);
        Warped_Moving_Mask_Hybrid1=zeros(d1,d2,d3);
        Warped_Moving_Mask_Hybrid5A1CA=zeros(d1,d2,d3);
        Warped_Moving_Mask_Hybrid5A1C=zeros(d1,d2,d3);
        Warped_Moving_Mask_Demons005=zeros(d1,d2,d3);
        Warped_Moving_Mask_Demons5=zeros(d1,d2,d3);

        for FN = 40:70
            Warped_Moving_Mask_Hybrid01(:,:,FN) = IHybrid01(:,:,FN);
            Warped_Moving_Mask_Hybrid2(:,:,FN) = IHybrid2(:,:,FN);
            Warped_Moving_Mask_Hybrid3(:,:,FN) = IHybrid3(:,:,FN);
            Warped_Moving_Mask_Hybrid5(:,:,FN) = IHybrid5(:,:,FN);
            Warped_Moving_Mask_Hybrid65(:,:,FN) = IHybrid65(:,:,FN);
            Warped_Moving_Mask_Hybrid75(:,:,FN) = IHybrid75(:,:,FN);
            Warped_Moving_Mask_Hybrid1(:,:,FN) = IHybrid1(:,:,FN);
            Warped_Moving_Mask_Hybrid5A1CA(:,:,FN) = IHybrid5A1CA(:,:,FN);
            Warped_Moving_Mask_Hybrid5A1C(:,:,FN) = IHybrid5A1C(:,:,FN);
            Warped_Moving_Mask_Demons005(:,:,FN) = IIODemons005(:,:,FN);
            Warped_Moving_Mask_Demons5(:,:,FN) = IIODemons5(:,:,FN);

            clear FN;            
        end        
        clear IHybrid01 IHybrid2 IHybrid3 IHybrid5 IHybrid65 IHybrid75 IHybrid1 IHybrid5A1CA IHybrid5A1C Hybrid5A2S IHybrid5A3S IIODemons005 IIODemons5;
        
        
        % DSC Calculation
        DSC_VH01(1,(MN-1)*2+SN) = DSC(Warped_Moving_Mask_Hybrid01,Fixed_Mask);
        DSC_VH2(1,(MN-1)*2+SN) = DSC(Warped_Moving_Mask_Hybrid2,Fixed_Mask);
        DSC_VH3(1,(MN-1)*2+SN) = DSC(Warped_Moving_Mask_Hybrid3,Fixed_Mask);
        DSC_VH5(1,(MN-1)*2+SN) = DSC(Warped_Moving_Mask_Hybrid5,Fixed_Mask);
        DSC_VH65(1,(MN-1)*2+SN) = DSC(Warped_Moving_Mask_Hybrid65,Fixed_Mask);
        DSC_VH75(1,(MN-1)*2+SN) = DSC(Warped_Moving_Mask_Hybrid75,Fixed_Mask);
        DSC_VH1(1,(MN-1)*2+SN) = DSC(Warped_Moving_Mask_Hybrid1,Fixed_Mask);
        DSC_VH5A1CA(1,(MN-1)*2+SN) = DSC(Warped_Moving_Mask_Hybrid5A1CA,Fixed_Mask);
        DSC_VH5A1C(1,(MN-1)*2+SN) = DSC(Warped_Moving_Mask_Hybrid5A1C,Fixed_Mask);
        DSC_VD005(1,(MN-1)*2+SN) = DSC(Warped_Moving_Mask_Demons005,Fixed_Mask);
        DSC_VD5(1,(MN-1)*2+SN) = DSC(Warped_Moving_Mask_Demons5,Fixed_Mask);
        
        %SN
        clear SN Fixed_Mask Warped_Moving_Mask_Hybrid01 Warped_Moving_Mask_Hybrid2 Warped_Moving_Mask_Hybrid3 Warped_Moving_Mask_Hybrid5 Warped_Moving_Mask_Hybrid65 Warped_Moving_Mask_Hybrid75 Warped_Moving_Mask_Hybrid1 Warped_Moving_Mask_Hybrid5A1CA Warped_Moving_Mask_Hybrid5A1C Warped_Moving_Mask_Demons005 Warped_Moving_Mask_Demons005;
    end % SN
    %MN
    clear MN;
end % MN
%Time=toc;
VH01=DSC_VH01;
VH2=DSC_VH2;
VH3=DSC_VH3;
VH5=DSC_VH5;
VH65=DSC_VH65;
VH75=DSC_VH75;
VH1=DSC_VH1;
VH5A1CA=DSC_VH5A1CA;
VH5A1C=DSC_VH5A1C;
VD005=DSC_VD005;
VD5=DSC_VD5;

clear  DSC_VH01 DSC_VH2 DSC_VH3 DSC_VH5 DSC_VH65 DSC_VH75 DSC_VH1 DSC_VH5A1CA DSC_VH5A1C DSC_VD005 DSC_VD5;

%mean(DSC_Vector)






