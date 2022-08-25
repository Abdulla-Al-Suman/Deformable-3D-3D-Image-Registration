clc;
close all
clear;
warning off all;
X=0:0.0001:127;
X(find(X>24)) = X(find(X>24))-48*round(X(find(X>24))/48);
A=BSpline_wavelet_S48(X);
a=0:0.0001:127;
plot(a,A,'LineWidth',2)
title('Dilated Periodic Spline Wavelet')
%xlabel('x')
%grid on
% ylabel('X')

 grid on
 axis([0 130 -0.2 1])
 set(gca,'FontName','Times New Roman','FontSize',15)
 set(gca,'XTick',0:20:130);
 set(gca,'YTick',-0.2:0.2:1);
 xlabel('x','FontName','Times New Roman','fontsize',23)
