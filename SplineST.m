function S = SplineST(X,jj,k)% X=1*L vector, SplineST= Spline Scalling translation, jj= scalling up , k= translation
%Fourth order B-Spline function
[~,L]=size(X);
S=zeros(1,L);
SSV=0;
i=0;
for x=X(1):X(2)-X(1):X(L)
    for j=0:4
        %a=(2^-6*x)-j; 
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


