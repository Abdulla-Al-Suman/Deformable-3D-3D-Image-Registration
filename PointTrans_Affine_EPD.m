function [xn,yn,zn] = PointTrans_Affine_EPD(m,x,y,z)

for i = 1:6
    xn{i} = m(1)*x{i}+m(2)*y{i}+m(3)*z{i}+m(4);
    yn{i} = m(5)*x{i}+m(6)*y{i}+m(7)*z{i}+m(8);
    zn{i} = m(9)*x{i}+m(10)*y{i}+m(11)*z{i}+m(12);
end

