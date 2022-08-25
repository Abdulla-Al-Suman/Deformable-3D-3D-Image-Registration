%function [Volumetric_DSCs_Hybrid65,Volumetric_DSCs_Hybrid1,Volumetric_DSCs_Demons] = Vol_LabelT_Warp_DSC_Cal_DT(F_muscle_data,M_muscle_data,d1,d2,d3,xOC,yOC,zOC,RP_WindowMI,RP_WindowEPD,RP_Affine19,RP_FFD5,RP_Affine087,RP_FFD4,RP_Affine044,RP_FFD3,RP_FFD2,RP_FFD1,DisH,DisD,dxdm5,dxdm4,dxdm3,dxdm2,dxdm1,NB5,NB4,NB3,NB2,NB1)
function [Volumetric_DSCs_Hybrid65,Volumetric_DSCs_Hybrid1,Volumetric_DSCs_Demons] = Vol_LabelT_Warp_DSC_Cal_DT(F_muscle_data,M_muscle_data,d1,d2,d3,xOC,yOC,zOC,RP_Affine19,RP_DPSW5_19,RP_FFD5,RP_FFD4,RP_FFD3,RP_FFD2,RP_FFD1,DisH,DisD,dxdmS96,dxdm5,dxdm4,dxdm3,dxdm2,dxdm1,BFD,NB5,NB4,NB3,NB2,NB1)
% transfering Labels using volume template process

% Creating vector for holding DSC values
DSC_V_Hybrid65 = zeros(1,6);
DSC_V_Hybrid1 = zeros(1,6);
DSC_V_Demons = zeros(1,6);


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
        I_ODemons = imwarp(Moving,DisD);
        I_ODemons = ImageThresholding(I_ODemons,0.02,d1,d2,d3);
        
        % START Window Search
        % data transfer for Window MI transformation
%         IW_MI = warp_Window_ST_3D(Moving,RP_WindowMI,xOC,yOC,zOC);
%         IW_MI65 = ImageThresholding(IW_MI,0.65,d1,d2,d3);
%         IW_MI1 = ImageThresholding(IW_MI,1,d1,d2,d3);
%         clear Moving;       
%       % data transfer for Window EPD transformation
%         IW_EPD65 = warp_Window_ST_3D(IW_MI65,RP_WindowEPD,xOC,yOC,zOC);
%         IW_EPD65 = ImageThresholding(IW_EPD65,0.65,d1,d2,d3);
%         IW_EPD1 = warp_Window_ST_3D(IW_MI1,RP_WindowEPD,xOC,yOC,zOC);
%         IW_EPD1 = ImageThresholding(IW_EPD1,1,d1,d2,d3);
%         clear IW_MI65 IW_MI1;
        % END Window Search
        
        % START 1st Section
        % data transfer for affine transformation
        IAffine19 = warp_AffineEPD_3D(Moving,RP_Affine19,xOC,yOC,zOC);
        IAffine19_65 = ImageThresholding(IAffine19,0.65,d1,d2,d3);
        IAffine19_1 = ImageThresholding(IAffine19,1,d1,d2,d3);        
        clear Moving IAffine19;        
        % data transfer for DPSW5 transformation
        IDPSW5_19_65 = warp_DPSWEPD_3D(RP_DPSW5_19,IAffine19_65,BFD,'linear',dxdmS96,xOC,yOC,zOC,d1,d2,d3);
        IDPSW5_19_65 = ImageThresholding(IDPSW5_19_65,0.65,d1,d2,d3);
        IDPSW5_19_1 = warp_DPSWEPD_3D(RP_DPSW5_19,IAffine19_1,BFD,'linear',dxdmS96,xOC,yOC,zOC,d1,d2,d3);
        IDPSW5_19_1 = ImageThresholding(IDPSW5_19_1,1,d1,d2,d3);
        clear IAffine19_65 IAffine19_1;
        % data transfer for DPSW4 transformation
%         IDPSW4_19 = warp_DPSWEPD_3D(RP_DPSW4_19,IDPSW5_19,BFD,'linear',dxdmS48,xOC,yOC,zOC,d1,d2,d3);
%         IDPSW4_19 = ImageThresholding(IDPSW4_19,0.65,d1,d2,d3);
%         clear IDPSW5_19;
         % data transfer for FFD5 transformation
        IFFD5_65=warp_FFDEPD_3D(RP_FFD5,IDPSW5_19_65,NB5,'linear',dxdm5,xOC,yOC,zOC,d1,d2,d3);
        IFFD5_65 = ImageThresholding(IFFD5_65,0.65,d1,d2,d3);
        IFFD5_1=warp_FFDEPD_3D(RP_FFD5,IDPSW5_19_1,NB5,'linear',dxdm5,xOC,yOC,zOC,d1,d2,d3);
        IFFD5_1 = ImageThresholding(IFFD5_1,1,d1,d2,d3);
        clear IDPSW5_19_65 IDPSW5_19_1;
        
        % START 2nd Section
        % data transfer for affine transformation
%         IAffine087_65 = warp_AffineEPD_3D(IFFD5_65,RP_Affine087,xOC,yOC,zOC);
%         IAffine087_65 = ImageThresholding(IAffine087_65,0.65,d1,d2,d3);
%         IAffine087_1 = warp_AffineEPD_3D(IFFD5_1,RP_Affine087,xOC,yOC,zOC);
%         IAffine087_1 = ImageThresholding(IAffine087_1,1,d1,d2,d3); 
%         clear IFFD5_65 IFFD5_1;        
%         % data transfer for DPSW5 transformation
%         IDPSW5_087_65 = warp_DPSWEPD_3D(RP_DPSW5_087,IAffine087_65,BFD,'linear',dxdmS96,xOC,yOC,zOC,d1,d2,d3);
%         IDPSW5_087_65 = ImageThresholding(IDPSW5_087_65,0.65,d1,d2,d3);
%         IDPSW5_087_1 = warp_DPSWEPD_3D(RP_DPSW5_087,IAffine087_1,BFD,'linear',dxdmS96,xOC,yOC,zOC,d1,d2,d3);
%         IDPSW5_087_1 = ImageThresholding(IDPSW5_087_1,1,d1,d2,d3);
%         clear IAffine087_65 IAffine087_1;
        % data transfer for DPSW4 transformation
%         IDPSW4_087_65 = warp_DPSWEPD_3D(RP_DPSW4_087,IDPSW5_087_65,BFD,'linear',dxdmS48,xOC,yOC,zOC,d1,d2,d3);
%         IDPSW4_087_65 = ImageThresholding(IDPSW4_087_65,0.65,d1,d2,d3);
%         IDPSW4_087_1 = warp_DPSWEPD_3D(RP_DPSW4_087,IDPSW5_087_1,BFD,'linear',dxdmS48,xOC,yOC,zOC,d1,d2,d3);
%         IDPSW4_087_1 = ImageThresholding(IDPSW4_087_1,1,d1,d2,d3);
%         clear IDPSW5_087_65 IDPSW5_087_1;
        % data transfer for FFD54 transformation
        IFFD4_65=warp_FFDEPD_3D(RP_FFD4,IFFD5_65,NB4,'linear',dxdm4,xOC,yOC,zOC,d1,d2,d3);
        IFFD4_65 = ImageThresholding(IFFD4_65,0.65,d1,d2,d3);
        IFFD4_1=warp_FFDEPD_3D(RP_FFD4,IFFD5_1,NB4,'linear',dxdm4,xOC,yOC,zOC,d1,d2,d3);
        IFFD4_1 = ImageThresholding(IFFD4_1,1,d1,d2,d3);
        clear IFFD5_65 IFFD5_1;
        
        % START 3rd Section
        % data transfer for affine transformation
%         IAffine044_65 = warp_AffineEPD_3D(IFFD4_65,RP_Affine044,xOC,yOC,zOC);
%         IAffine044_65 = ImageThresholding(IAffine044_65,0.65,d1,d2,d3); 
%         IAffine044_1 = warp_AffineEPD_3D(IFFD4_1,RP_Affine044,xOC,yOC,zOC);
%         IAffine044_1 = ImageThresholding(IAffine044_1,1,d1,d2,d3);
%         clear IFFD4_65 IFFD4_1;        
%         % data transfer for DPSW5 transformation
%         IDPSW5_044_65 = warp_DPSWEPD_3D(RP_DPSW5_044,IAffine044_65,BFD,'linear',dxdmS96,xOC,yOC,zOC,d1,d2,d3);
%         IDPSW5_044_65 = ImageThresholding(IDPSW5_044_65,0.65,d1,d2,d3);
%         IDPSW5_044_1 = warp_DPSWEPD_3D(RP_DPSW5_044,IAffine044_1,BFD,'linear',dxdmS96,xOC,yOC,zOC,d1,d2,d3);
%         IDPSW5_044_1 = ImageThresholding(IDPSW5_044_1,1,d1,d2,d3);
%         clear IAffine044_65 IAffine044_1;
        % data transfer for DPSW4 transformation
%         IDPSW4_044_65 = warp_DPSWEPD_3D(RP_DPSW4_044,IDPSW5_044_65,BFD,'linear',dxdmS48,xOC,yOC,zOC,d1,d2,d3);
%         IDPSW4_044_65 = ImageThresholding(IDPSW4_044_65,0.65,d1,d2,d3);
%         IDPSW4_044_1 = warp_DPSWEPD_3D(RP_DPSW4_044,IDPSW5_044_1,BFD,'linear',dxdmS48,xOC,yOC,zOC,d1,d2,d3);
%         IDPSW4_044_1 = ImageThresholding(IDPSW4_044_1,1,d1,d2,d3);
%         clear IDPSW5_044_65 IDPSW5_044_1;
        % data transfer for FFD543 transformation
        IFFD3_65=warp_FFDEPD_3D(RP_FFD3,IFFD4_65,NB3,'linear',dxdm3,xOC,yOC,zOC,d1,d2,d3);
        IFFD3_65 = ImageThresholding(IFFD3_65,0.65,d1,d2,d3);
        IFFD3_1=warp_FFDEPD_3D(RP_FFD3,IFFD4_1,NB3,'linear',dxdm3,xOC,yOC,zOC,d1,d2,d3);
        IFFD3_1 = ImageThresholding(IFFD3_1,1,d1,d2,d3);
        clear IFFD4_65 IFFD4_1;
        
        % START 4th Section
        % data transfer for affine transformation
%         IAffine012 = warp_AffineEPD_3D(IFFD3,RP_Affine012,xOC,yOC,zOC);
%         IAffine012 = ImageThresholding(IAffine012,0.65,d1,d2,d3); 
%         clear IFFD3;        
%         % data transfer for DPSW5 transformation
%         IDPSW5_012 = warp_DPSWEPD_3D(R_P_DPSW5_012,IAffine012,BFD,'linear',dxdmS96,xOC,yOC,zOC,d1,d2,d3);
%         IDPSW5_012 = ImageThresholding(IDPSW5_012,0.65,d1,d2,d3);
%         clear IAffine012;
%         % data transfer for DPSW4 transformation
%         IDPSW4_012 = warp_DPSWEPD_3D(R_P_DPSW4_012,IDPSW5_012,BFD,'linear',dxdmS48,xOC,yOC,zOC,d1,d2,d3);
%         IDPSW4_012 = ImageThresholding(IDPSW4_012,0.65,d1,d2,d3);
%         clear IDPSW5_012;
        % data transfer for FFD5432 transformation
        IFFD2_65=warp_FFDEPD_3D(RP_FFD2,IFFD3_65,NB2,'linear',dxdm2,xOC,yOC,zOC,d1,d2,d3);
        IFFD2_65 = ImageThresholding(IFFD2_65,0.65,d1,d2,d3);
        IFFD2_1=warp_FFDEPD_3D(RP_FFD2,IFFD3_1,NB2,'linear',dxdm2,xOC,yOC,zOC,d1,d2,d3);
        IFFD2_1 = ImageThresholding(IFFD2_1,1,d1,d2,d3);
        clear IFFD3_65 IFFD3_1;
        
        % START 5th Section
        % data transfer for affine transformation
%         IAffine011 = warp_AffineEPD_3D(IFFD2,R_P_Affine011,xOC,yOC,zOC);
%         IAffine011 = ImageThresholding(IAffine011,0.65,d1,d2,d3); 
%         clear IFFD2;        
%         % data transfer for DPSW5 transformation
%         IDPSW5_011 = warp_DPSWEPD_3D(R_P_DPSW5_011,IAffine011,BFD,'linear',dxdmS96,xOC,yOC,zOC,d1,d2,d3);
%         IDPSW5_011 = ImageThresholding(IDPSW5_011,0.65,d1,d2,d3);
%         clear IAffine011;
%         % data transfer for DPSW4 transformation
%         IDPSW4_011 = warp_DPSWEPD_3D(R_P_DPSW4_011,IDPSW5_011,BFD,'linear',dxdmS48,xOC,yOC,zOC,d1,d2,d3);
%         IDPSW4_011 = ImageThresholding(IDPSW4_011,0.65,d1,d2,d3);
%         clear IDPSW5_011;
        % data transfer for FFD54321 transformation
        IFFD1_65=warp_FFDEPD_3D(RP_FFD1,IFFD2_65,NB1,'linear',dxdm1,xOC,yOC,zOC,d1,d2,d3);
        IFFD1_65 = ImageThresholding(IFFD1_65,0.65,d1,d2,d3);
        IFFD1_1=warp_FFDEPD_3D(RP_FFD1,IFFD2_1,NB1,'linear',dxdm1,xOC,yOC,zOC,d1,d2,d3);
        IFFD1_1 = ImageThresholding(IFFD1_1,1,d1,d2,d3);
        clear IFFD2_65 IFFD2_1;
        
        % data transfer for Hybrid-Demon
        IHybrid65 = imwarp(IFFD1_65,DisH);
        IHybrid65 = ImageThresholding(IHybrid65,.02,d1,d2,d3);
        IHybrid1 = imwarp(IFFD1_1,DisH);
        IHybrid1 = ImageThresholding(IHybrid1,.02,d1,d2,d3);
        clear IFFD1_65 IFFD1_1;

        
        % making ROI warped moving Image
        Warped_Moving_Mask_Hybrid65=zeros(d1,d2,d3);
        Warped_Moving_Mask_Hybrid1=zeros(d1,d2,d3);
        Warped_Moving_Mask_Demons=zeros(d1,d2,d3);

        for FN = 40:70
            Warped_Moving_Mask_Hybrid65(:,:,FN) = IHybrid65(:,:,FN);
            Warped_Moving_Mask_Hybrid1(:,:,FN) = IHybrid1(:,:,FN);
            Warped_Moving_Mask_Demons(:,:,FN) = I_ODemons(:,:,FN);
            clear FN;            
        end        
        clear IHybrid65 IHybrid1 I_ODemons;
        
        
        % DSC Calculation
        %DSC_V_Hybrid(1,(MN-1)*2+SN) = dice(Warped_Moving_Mask_Hybrid,Fixed_Mask);
        DSC_V_Hybrid65(1,(MN-1)*2+SN) = DSC(Warped_Moving_Mask_Hybrid65,Fixed_Mask);
        DSC_V_Hybrid1(1,(MN-1)*2+SN) = DSC(Warped_Moving_Mask_Hybrid1,Fixed_Mask);
        DSC_V_Demons(1,(MN-1)*2+SN) = DSC(Warped_Moving_Mask_Demons,Fixed_Mask);
        
        %SN
        clear SN Fixed_Mask Warped_Moving_Mask_Hybrid65 Warped_Moving_Mask_Hybrid1 Warped_Moving_Mask_Demons;
    end % SN
    %MN
    clear MN;
end % MN
%Time=toc;
Volumetric_DSCs_Hybrid65=DSC_V_Hybrid65;
Volumetric_DSCs_Hybrid1=DSC_V_Hybrid1;
Volumetric_DSCs_Demons=DSC_V_Demons;

clear DSC_V_Hybrid65 DSC_V_Hybrid1 DSC_V_Demons;

%mean(DSC_Vector)






