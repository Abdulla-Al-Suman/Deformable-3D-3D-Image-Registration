
function [Smin] = update_m_DPSWEPD(Smin,I,edge_list_R,iter,d1,d2,d3,s,dxdm,LT,UT,n)


global sm mi m
EI = zeros(d1,d2,d3);
for i=1:d3
EI(:,:,i) = edge(I(:,:,i),'canny',[LT UT],1.5); % Calculation of Canny Edge detection of moving image volume
end
DI = ChamferDis3dinfv2(uint32(EI),inf);
[Gx,Gy,Gz] = gradient(DI);

S = CalcEPD3D(edge_list_R,DI);
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

% dPdm calculation START**********
%dPdm = calc_dPdm(Gx,Gy,Gz,s,dxdm,d1,d2,d3);
dPdm  = zeros(d1,d2,d3,3*s^3);

for k=1:s^3
    
    dPdm(:,:,:,k) = Gx.*squeeze(dxdm(:,:,:,k));
    dPdm(:,:,:,k+s^3) = Gy.*squeeze(dxdm(:,:,:,k));
    dPdm(:,:,:,k+(2*s^3)) = Gz.*squeeze(dxdm(:,:,:,k));
end
clear Gx Gy Gz k;
% dPdm calculation END**********

%***** m update  calculation START***************
%m = update_mm(dPdm,m,edge_list_R,s);
J = zeros(3*s^3,n);
b = zeros(1,3*s^3);

for i = 1:3*s^3
    B = 0;
    for j = 1:n
        J(i,j) = squeeze(dPdm(edge_list_R(j,2),edge_list_R(j,1),edge_list_R(j,3),i));
        B = B+squeeze(dPdm(edge_list_R(j,2),edge_list_R(j,1),edge_list_R(j,3),i));
    end
     b(i) = B;
     clear B;
end

H = J*J';
dm = -inv(H)*b';
m  = m+dm';
%***** m update  calculation END***************
sm(iter,:)=m;
clear dPdm i j b B H J dm;

% function m = update_mm(dPdm,m,edge_list,s)
% 
% [n,~,~] = size(edge_list);
% J = zeros(3*s^3,n);
% b = zeros(1,3*s^3);
% 
% for i = 1:3*s^3
%     B = 0;
%     for j = 1:n
%         J(i,j) = squeeze(dPdm(edge_list(j,2),edge_list(j,1),edge_list(j,3),i));
%         B = B+squeeze(dPdm(edge_list(j,2),edge_list(j,1),edge_list(j,3),i));
%     end
%      b(i) = B;
% end
% 
% H = J*J';
% dm = -inv(H)*b';
% m  = m+dm';

% function dPdm = calc_dPdm(Gx,Gy,Gz,s,dxdm,d1,d2,d3)
% 
% dPdm  = zeros(d1,d2,d3,3*s^3);
% 
% for k=1:s^3
%     
%     dPdm(:,:,:,k) = Gx.*squeeze(dxdm(:,:,:,k));
%     dPdm(:,:,:,k+s^3) = Gy.*squeeze(dxdm(:,:,:,k));
%     dPdm(:,:,:,k+(2*s^3)) = Gz.*squeeze(dxdm(:,:,:,k));
% end
% clear Gx Gy Gz k;


