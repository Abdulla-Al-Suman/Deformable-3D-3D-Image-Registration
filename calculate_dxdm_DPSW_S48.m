function dxdmS48 = calculate_dxdm_DPSW_S48(s,d1,d2,d3,x,y,z)

dxdmS48 = zeros(d1,d2,d3,s^3);
k = 1;
for u = 0:s-1
    for v = 0:s-1
        for w=0:s-1
            XW=x*u;
            XW(find(XW>24)) = XW(find(XW>24))-48*round(XW(find(XW>24))/48); 
            YW=y*v;
            YW(find(YW>24)) = YW(find(YW>24))-48*round(YW(find(YW>24))/48);
            ZW=z*w;
            ZW(find(ZW>24)) = ZW(find(ZW>24))-48*round(ZW(find(ZW>24))/48);
            dxdmS48(:,:,:,k) = BSpline_wavelet_S48(XW).*BSpline_wavelet_S48(YW).*BSpline_wavelet_S48(ZW);
            k = k+1;
            clear  XW YW ZW;
        end
    end
end
clear k u v w;