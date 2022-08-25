function [fxn,fyn,fzn] = PointTrans_Demons(D,X,Y,Z,d1,d2,d3)

for i = 1:6
    X_Dis=D(:,:,:,1);
    Y_Dis=D(:,:,:,2);
    Z_Dis=D(:,:,:,2);
    fxn{i} = X{i}+X_Dis(sub2ind(size(X_Dis),round(Y{i}+(d1/2+0.5)),round(X{i}+(d2/2+0.5)),round(Z{i}+(d3/2+0.5))));%transfer center2Corner for matrix access
    fyn{i} = Y{i}+Y_Dis(sub2ind(size(Y_Dis),round(Y{i}+(d1/2+0.5)),round(X{i}+(d2/2+0.5)),round(Z{i}+(d3/2+0.5))));
    fzn{i} = Z{i}+Z_Dis(sub2ind(size(Z_Dis),round(Y{i}+(d1/2+0.5)),round(X{i}+(d2/2+0.5)),round(Z{i}+(d3/2+0.5))));
    clear X_Dis Y_Dis Z_Dis;
end
