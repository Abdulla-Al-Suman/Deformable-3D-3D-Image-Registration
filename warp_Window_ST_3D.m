function V = warp_Window_ST_3D(V0,m,x,y,z)
%global  m
% [d1,d2,d3]=size(V0);
% [x,y,z] = meshgrid((0:d2-1)-d2/2,(0:d1-1)-d1/2,(0:d3-1)-d3/2);% x y z generate again

%M = [m(1) m(2) m(3) m(4); m(5) m(6) m(7)  m(8);m(9)  m(10) m(11) m(12); 0 0 0 1];
%M = makeM_affine_3D(m(1),m(2),m(3),m(4),m(5),m(6),m(7),m(8),m(9),m(10),m(11),m(12),m(13),m(14),m(15));  
xn = m(1)*x+m(2);
yn = m(3)*y+m(4);
zn = m(5)*z+m(6);

% xn = M(1,1)*x+M(1,2)*y+M(1,3)*z+M(1,4);
% yn = M(2,1)*x+M(2,2)*y+M(2,3)*z+M(2,4);
% zn = M(3,1)*x+M(3,2)*y+M(3,3)*z+M(3,4);
%V=interp3(x,y,z,V0,xn,yn,zn,method,0); 
%V=interp3_mex(x,y,z,V0,xn,yn,zn);
V = interp3(x,y,z,V0,xn,yn,zn,'linear',0);
%V = interp3(x,y,z,V0,xn,yn,zn,'spline');
clear xn yn zn;




