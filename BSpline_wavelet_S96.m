function DSW = BSpline_wavelet_S96(x)
%function  BSpline_wavelet(x)
% clc
% clear
% *********Start Scalling function support[0 4] ***************
MX=4;%MX=maximum X value
SS=0.0001; % SS= step size should be .1,.2,.5 i.e dividide 1 into equal parts 
X=0:SS:MX;
[~,L1]=size(X);
SF=Scalling_Function(X,SS,MX); % Scalling function support[0 4]
% figure(1)
% plot(X,SF,'r')
% *********END Scalling function support[0 4] ***************

% *********Start Spline wavelet generation [0 4]***************
SF1 = padarray(SF,[0 2*L1-L1],'post');

SF2=padarray(SF,[0 (1/SS)],'pre'); %1/SS= no of step in 1, right time shifting operation
[~,L2]=size(SF2);
SF2 = padarray(SF2,[0 2*L1-L2],'post');

SF3=padarray(SF,[0 (2/SS)],'pre');
[~,L3]=size(SF3);
SF3 = padarray(SF3,[0 2*L1-L3],'post');

SW=-(3/7)*SF1(round((2*X*(1/SS))+1))+(12/7)*SF2(round((2*X*(1/SS))+1))-(3/7)*SF3(round((2*X*(1/SS))+1));
%DSW=-(3/7)*SF1(round((2*x*(1/SS))+1))+(12/7)*SF2(round((2*x*(1/SS))+1))-(3/7)*SF3(round((2*x*(1/SS))+1));

% figure(2)
% plot(X,SW,'r')

X192=0:SS:192;%X1p5= 1 point 5 support
SW192=SW(round((X192/(SS*64))+1));
%DSW=SW192(round(((x+96)*(1/SS))+1));
X96=0:SS:96; % Have to change
SW96=SW192(round(((X96*2)/SS)+1)); % scalling up means divided. scalling down meand multiplication
DSW=SW96(round(((x+48)*(1/SS))+1));

% x=0:.0001:192;
% plot(x,BSpline_wavelet(x))

%**************************** END Spline wavelet generation [0 4] **********

% ******** Mapping Scalling function [0 4] to [0 256] ************

%     SF256=SF(round((x/(SS*64))+1)); % Scalling function support[0 256]

% ******** End Mapping Scalling function [0 4] to [0 256] ************

% % ******** Mapping Wavelet function [0 4] to [0 256] ************
% X192=0:SS:192;
% SW192=SW(round((X192/(SS*64))+1)); % 192 support per wavelet
% % ******** End Mapping Wavelet function [0 4] to [0 256] ************
% 
% % ******** Generation 4th dilated  Wavelet function,   ************
% X12=0:SS:12;
% SW12=SW192(round((128*X12*(1/SS))+1)); %12 support per wavelet
% % ******** End 4th dilated  Wavelet function ************
% 
% % ******** Generation 6 right shift of 4th dilated  Wavelet function,   ************
% %SX12=x-12;
% %SX12=-6:SS:6; % SX12 = shifted X12, Period or support=12
% % if x<=6
% %     SW12(round(((x+6)*(1/SS))+1))
%     DSW=SW12(round(((x+6)*(1/SS))+1)); %12 support per wavelet
% %end
% % elseif x>6
%      %DSW=SW12(round(((x-6)*(1/SS))+1));
% % end
% ******** End 6 right shift of 4th dilated  Wavelet function ************

% % ******** Generation horizontal stacking of Wavelets,   ************
% DSW=SSW12(round(((x+6)*(1/SS))+1)); %12 support per wavelet
% % ******** End horizontal stacking of Wavelets, ************

% Generation wavelet basis function by translating over -127 to 256
% XXX=-127:.01:257;k
% [~,L4]=size(XXX);
% FSW=zeros(1,L4);%FSW=Final spline wavelet
% for j=1:128
%     if j==1
%         FSW(1:301) = SW;
%     else
%         FSW((j-1)*300+2:j*300+1) = SW(2:301);
%     end
% end
% figure(3)
% plot(XXX,FSW)







