function [R_Moving_V_WindowMI,R_Parameters_WindowMI] = Registration_Window_MI(Fixed_Volume,Moving_Volume,x,y,z)

% initialize EPD similarity measure
Smax = 0;
Fixed_Volume_W = Fixed_Volume - min(Fixed_Volume(:));
Fixed_Volume_W  = round(255.4*Fixed_Volume_W /(max(Fixed_Volume_W (:)))); 
for tx = -8:2:8
    for ty = -14:2:14 %-w-20:4:w+20        
        R_Parameters_WindowMI=[1.0000 tx 1.0000 ty 1.0000 0.0000];
        R_Moving_V_WindowMI = warp_Window_ST_3D(Moving_Volume,R_Parameters_WindowMI,x,y,z);
        R_Moving_V_WindowMI = R_Moving_V_WindowMI - min(R_Moving_V_WindowMI(:));
        R_Moving_V_WindowMI = round(255.4*R_Moving_V_WindowMI/(max(R_Moving_V_WindowMI(:))));
        
        
        [~,~,S] = ent(Fixed_Volume_W,R_Moving_V_WindowMI);
        if S > Smax
            Smax = S;
            txm = tx;
            tym = ty;
        end
        clear ty S R_Parameters_WindowMI R_Moving_V_WindowMI;
    end
    clear tx;
end
R_Parameters_WindowMI=[1.0000 txm 1.0000 tym 1.0000 0.0000];
R_Moving_V_WindowMI = warp_Window_ST_3D(Moving_Volume,R_Parameters_WindowMI,x,y,z);
clear Smax txm tym;
%clear Smax tym;


