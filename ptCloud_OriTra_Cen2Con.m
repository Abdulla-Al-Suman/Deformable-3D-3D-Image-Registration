function [Xt,Yt,Zt]=ptCloud_OriTra_Cen2Con(X,Y,Z,d1,d2,d3)

for i =1:6
Xt{i}=X{i}+(d2/2+0.5);Yt{i}=Y{i}+(d1/2+0.5);Zt{i}=Z{i}+(d3/2);
end

