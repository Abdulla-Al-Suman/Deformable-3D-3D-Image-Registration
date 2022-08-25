function DSW = BSpline_wavelet_S48(x)
% *********Start Scalling function support[0 4] ***************
MX=4;%MX=maximum X value
SS=0.0001; % SS= step size should be .1,.2,.5 i.e divide 1 into equal parts 
X=0:SS:MX;
[~,L1]=size(X);
SF=Scalling_Function(X,SS,MX); % Scalling function support[0 4]
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
% *********END Spline wavelet generation [0 4]***************
% 4th dilated with left shifted wavelet generation
X48=0:SS:48;
SW48=SW(round((X48/(SS*16))+1));
DSW=SW48(round(((x+24)*(1/SS))+1));

