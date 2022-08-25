function V = warp_DPSWEPD_3D(M,V0,s,method,dxdm,x,y,z,d1,d2,d3)


xn = zeros(d1,d2,d3);
yn = zeros(d1,d2,d3);
zn = zeros(d1,d2,d3);

for k = 1:s^3
    xn = xn+M(k)*squeeze(dxdm(:,:,:,k));
    yn = yn+M(k+s^3)*squeeze(dxdm(:,:,:,k));
    zn = zn+M(k+(2*s^3))*squeeze(dxdm(:,:,:,k));
end

fxn = x+xn;
fyn = y+yn;
fzn = z+zn;

%  V = interp3(x,y,z,V0,xn,yn,zn,'linear',0);
V = interp3(x,y,z,V0,fxn,fyn,fzn,method,0);
clear xn yn zn fxn fyn fzn;

% figure(2), imagesc([V(:,:,5)  V0(:,:,5)]) 
