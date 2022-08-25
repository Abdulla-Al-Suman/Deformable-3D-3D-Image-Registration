function dxdmS96 = calculate_dxdm_DPSW_S96(s,d1,d2,d3,x,y,z)

dxdmS96 = zeros(d1,d2,d3,s^3);
k = 1;
for u = 0:s-1
    for v = 0:s-1
        for w=0:s-1
            XW=x*u;
            XW(find(XW>48)) = XW(find(XW>48))-96*round(XW(find(XW>48))/96); 
            YW=y*v;
            YW(find(YW>48)) = YW(find(YW>48))-96*round(YW(find(YW>48))/96);
            ZW=z*w;
            ZW(find(ZW>48)) = ZW(find(ZW>48))-96*round(ZW(find(ZW>48))/96);
            dxdmS96(:,:,:,k) = BSpline_wavelet_S96(XW).*BSpline_wavelet_S96(YW).*BSpline_wavelet_S96(ZW);
            k = k+1;
            clear  XW YW ZW;
        end
    end
end
clear k u v w;