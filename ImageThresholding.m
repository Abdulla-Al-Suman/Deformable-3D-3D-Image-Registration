function II = ImageThresholding(III,TH,d1,d2,d3)

II=zeros(d1,d2,d3);
for k=1:d1
    for i=1:d2
        for j=1:d3
            if III(i,j,k)>=TH
                II(i,j,k)=1;
%             else
%                 II(i,j,k)=0;
            end
            
        end
    end
end
clear  i j k;