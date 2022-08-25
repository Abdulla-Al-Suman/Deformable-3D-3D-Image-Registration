function dxdm = calculate_dxdm_FFD(d1,d2,d3,NB,jj)
dxdm = zeros(d1,d2,d3,NB);
X=1:d2;
ii=0;
for k=-((2^-jj)*d1)/2:((2^-jj)*d1)/2  % control points, translation index
    TensorP2 = zeros(d1,d2,d3);
    ii=ii+1;
    Spline = SplineST(X,jj,k);
    KronProduct1 = kron(Spline,Spline');
    for i=1:d3
        OrderPairs2=KronProduct1*Spline(i);
        TensorP2(:,:,i)=OrderPairs2;
        clear OrderPairs2;
    end
    dxdm(:,:,:,ii)=2^-(jj)*TensorP2;
    clear KronProduct1 i TensorP2 Spline;
end
clear X ii;

function S = SplineST(X,jj,k)% X=1*L vector, SplineST= Spline Scalling translation, jj= scalling up , k= translation
%Fourth order B-Spline function
[~,L]=size(X);
S=zeros(1,L);
SSV=0;
i=0;
for x=X(1):X(2)-X(1):X(L)
    for j=0:4 
        a=(2^-jj*x)-j;  % 2^jj , jj= resolution level, 6,5,4, 
        if a<0
            SSV=SSV+0;
        else
        SSV=SSV+nchoosek(4,j)*(-1)^j*(a)^3;
        end
    end
    i=i+1;
    ii=i-(k*2^jj); %interval=2^jj
    if ii<1 || ii>L
         SSV=0;
    else
    S(1,ii)=(1/6)*SSV;
     SSV=0;
    end    
end