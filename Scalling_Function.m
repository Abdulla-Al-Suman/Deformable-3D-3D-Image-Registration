function Sf = Scalling_Function(X,SS,MX)% X=1*L vector
%Fourth order B-Spline function
[~,L]=size(X);
Sf=zeros(1,L);
SSV=0;
i=0;
for x=0:SS:MX
    for j=0:4
        a=x-j;
        if a<0
            SSV=SSV+0;
        else
        SSV=SSV+nchoosek(4,j)*(-1)^j*(a)^3;
        end
    end
    i=i+1;
    Sf(1,i)=(1/6)*SSV;
    SSV=0;
end
