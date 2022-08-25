
function [Smin] = update_m_AffineEPD(I,Smin,edge_list,iter,d1,d2,d3,n,x,y,z,LT,UT)

global sm mi m
%[d1,d2,d3]=size(I);
%[x,y,z] = meshgrid((0:d1-1)-d1/2,(0:d2-1)-d2/2,(0:d3-1)-d3/2);

EI = zeros(d1,d2,d3);
for i=1:d3
EI(:,:,i) = edge(I(:,:,i),'canny',[LT UT],1.5); % Calculation of Canny Edge detection of moving image volume
end

% tic
 DI = ChamferDis3dinfv2(uint32(EI),inf);
% toc
%pause
[Gx,Gy,Gz] = gradient(DI);
S = CalcEPD3D(edge_list,DI);
%sepd(1,iter)=S;
if S < Smin
	   Smin = S; % update minimum EPD value if EPD value is less than previous minimum EPD
       
       if iter==1
           mi=0; %mi= minimum iteration
       else
           mi=iter-1;
       end    
end
clear S DI EI i;
%*****dPdm calculation START***************
dPdm = zeros(12,d1,d2,d3);% Find values for dI'/dm1, dI'/dm2, dI'/dm3 and  dI'/dm4

dPdm(1,:,:,:) = x.*Gx;
dPdm(2,:,:,:) = y.*Gx;
dPdm(3,:,:,:) = z.*Gx;
dPdm(4,:,:,:) = Gx;
% Find values for dI'/dm5, dI'/dm6, dI'/dm7 and dI'/dm8
dPdm(5,:,:,:) = x.*Gy;
dPdm(6,:,:,:) = y.*Gy;
dPdm(7,:,:,:) = z.*Gy;
dPdm(8,:,:,:) = Gy;
% Find values for dI'/dm9, dI'/dm10, dI'/dm11 and dI'/dm12
dPdm(9,:,:,:) = x.*Gz;
dPdm(10,:,:,:) = y.*Gz;
dPdm(11,:,:,:) = z.*Gz;
dPdm(12,:,:,:) = Gz;
clear Gx Gy Gz;
%*****dPdm calculation END***************
%***** m update  calculation START***************
J = zeros(12,n);
b = zeros(1,12);

for i = 1:12
    B = 0;
    for j = 1:n
        J(i,j) = squeeze(dPdm(i,edge_list(j,2),edge_list(j,1),edge_list(j,3)));
        B = B+squeeze(dPdm(i,edge_list(j,2),edge_list(j,1),edge_list(j,3)));
    end
     b(i) = B;
end

H = J*J';
dm = -inv(H)*b';
% dm = dm';
m  = m+dm';
%***** m update  calculation END***************

sm(iter,:)=m;
clear dPdm H dm i j J b B;



