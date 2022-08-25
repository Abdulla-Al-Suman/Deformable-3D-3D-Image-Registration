function [fxn,fyn,fzn] = PointTrans_DPSWEPD(M,X,Y,Z,dxdm,d1,d2,d3)

s = 2;

for i = 1:6
    xn = zeros(1,length(X{i}));
    yn = zeros(1,length(Y{i}));
    zn = zeros(1,length(Z{i}));
    
       
        for k = 1:s^3
            
            PS_BasisFunction=squeeze(dxdm(:,:,:,k));
            E_PS_BasisFunction=PS_BasisFunction(sub2ind(size(PS_BasisFunction),round(Y{i}+(d1/2+0.5)),round(X{i}+(d2/2+0.5)),round(Z{i}+d3/2)));%transfer center2Corner for matrix access
            xn = xn+M(k)*E_PS_BasisFunction;
            yn = yn+M(k+s^3)*E_PS_BasisFunction;
            zn = zn+M(k+(2*s^3))*E_PS_BasisFunction;
            clear PS_BasisFunction E_PS_BasisFunction;
        end
    
    fxn{i} = X{i}+xn;
    fyn{i} = Y{i}+yn;
    fzn{i} = Z{i}+zn;
    clear xn yn zn;
end

