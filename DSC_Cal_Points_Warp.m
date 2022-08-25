function [Volumetric_DSCs_Hybrid,Volumetric_DSCs_Demons] = DSC_Cal_Points_Warp(Xm,Ym,Zm,Xf_TH,Yf_TH,Zf_TH,Xf_TD,Yf_TD,Zf_TD,d1,d2,d3)
% transfering Fixed to Moving using 

% Creating vector for holding DSC values
DSC_V_Hybrid = zeros(1,6);
DSC_V_Demons = zeros(1,6);


muscle = [3 16 22]; % muscles numbers
for MN = 1:length(muscle) % MN=muscle number
    for SN = 1:2 %SNA=1(Left)= side number,
        
        Moving_Vol = ptCloud2vol(Xm,Ym,Zm,(MN-1)*2+SN);
        W_Fixed_Hybrid_Vol = ptCloud2vol(Xf_TH,Yf_TH,Zf_TH,(MN-1)*2+SN);
        W_Fixed_Demons_Vol = ptCloud2vol(Xf_TD,Yf_TD,Zf_TD,(MN-1)*2+SN);
        
        
        % making ROI warped Fixed mask
        Warped_Fixed_Mask_Hybrid=zeros(d1,d2,d3);
        Warped_Fixed_Mask_Demons=zeros(d1,d2,d3);
        Moving_Mask=zeros(d1,d2,d3);

        for FN = 40:70
            Warped_Fixed_Mask_Hybrid(:,:,FN) = W_Fixed_Hybrid_Vol(:,:,FN);
            Warped_Fixed_Mask_Demons(:,:,FN) = W_Fixed_Demons_Vol(:,:,FN);
            Moving_Mask(:,:,FN) = Moving_Vol(:,:,FN);
            clear FN;            
        end        
        clear W_Fixed_Hybrid_Vol W_Fixed_Demons_Vol Moving_Vol;
        
        
        % DSC Calculation
        %DSC_V_Hybrid(1,(MN-1)*2+SN) = dice(Warped_Moving_Mask_Hybrid,Fixed_Mask);
        DSC_V_Hybrid(1,(MN-1)*2+SN) = DSC(Warped_Fixed_Mask_Hybrid,Moving_Mask);
        DSC_V_Demons(1,(MN-1)*2+SN) = DSC(Warped_Fixed_Mask_Demons,Moving_Mask);
        
        %SN
        clear SN Moving_Mask Warped_Fixed_Mask_Hybrid Warped_Fixed_Mask_Demons;
    end % SN
    %MN
    clear MN;
end % MN
%Time=toc;
Volumetric_DSCs_Hybrid=DSC_V_Hybrid;
Volumetric_DSCs_Demons=DSC_V_Demons;
clear DSC_V_Hybrid DSC_V_Demons ;






