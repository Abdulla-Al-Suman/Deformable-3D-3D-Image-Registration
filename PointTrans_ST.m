function [xn,yn,zn] = PointTrans_ST(m,x,y,z)

for i = 1:6
    xn{i} = m(1)*x{i}+m(2);
    yn{i} = m(3)*y{i}+m(4);
    zn{i} = m(5)*z{i}+m(6);
end

