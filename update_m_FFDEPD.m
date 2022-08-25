
function [Smin] = update_m_FFDEPD(Smin,I,edge_list_R,iter,d1,d2,d3,NB,dxdm,n,LT,UT)

global sm mi m
EI = zeros(d1,d2,d3);
for i=1:d3
EI(:,:,i) = edge(I(:,:,i),'canny',[LT UT],1.5); % Calculation of Canny Edge detection of moving image volume
end
%DI = ChamferDis3d(uint32(EI));
DI = ChamferDis3dinfv2(uint32(EI),inf);
%figure(1), imagesc(DI(:,:,2)), colormap(gray)
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
%*****dPdm calculation START***************
dPdm  = zeros(d1,d2,d3,3*NB);
for k=1:NB   
    %dPdm(:,:,:,k) = single(Gx.*squeeze(dxdm(:,:,:,k)));
    dPdm(:,:,:,k) = Gx.*squeeze(dxdm(:,:,:,k));
    dPdm(:,:,:,k+NB) = Gy.*squeeze(dxdm(:,:,:,k));
    dPdm(:,:,:,k+(2*NB)) = Gz.*squeeze(dxdm(:,:,:,k));
end
clear Gx Gy Gz k;

%*****dPdm calculation END***************
%***** m update  calculation START***************
%J = zeros(3*NB,n);
b = zeros(1,3*NB);
for i = 1:3*NB
    B = 0;
    for j = 1:n
        %J(i,j) = squeeze(dPdm(edge_list_R(j,2),edge_list_R(j,1),edge_list_R(j,3),i));
        B = B+squeeze(dPdm(edge_list_R(j,2),edge_list_R(j,1),edge_list_R(j,3),i));
    end
     b(i) = B;
end

%H = J*J';
%dm = -inv(H)*b';

m  = m-b;
%***** m update  calculation END***************
sm(iter,:)=m;
clear dPdm i j b B;


