function [Xn,Yn,Zn] = get_ptCloud(muscle_data)

muscle = [3 16 22]; % muscles numbers
n = 1;
for MN = 1:length(muscle) % MN=muscle number
    for SN = 1:2 %SNA=1(Left)= side number,
        
        % START Point Cloud Generation
        X=zeros(1,1);Y=zeros(1,1);Z=zeros(1,1);
        for FN =1:128  % 1:O for full data, FN= frame number, for making banary volume
            x1 = muscle_data{muscle(MN)}{SN}{FN}.x;  % muscle selection, annotation software save x as 1 to N
            y1 = muscle_data{muscle(MN)}{SN}{FN}.y;  % 1= left, 2=right in the original annotation software
            L=length(x1);
            if L>=1
                %RESAMPLING
                x1(end+1) = x1(1);y1(end+1) = y1(1); % puting first element at end position, 
                
                x = 1:length(x1);v = x1;xq = 1:0.001:length(x1);x1=[];
                x1 = interp1(x,v,xq,'spline');
                x = 1:length(y1);v = y1;xq = 1:0.001:length(y1);y1=[];
                y1 = interp1(x,v,xq,'spline');
                L3 = length(x1);
                z1=FN*ones(1,L3);
                X = horzcat(X,x1);Y = horzcat(Y,y1);Z = horzcat(Z,z1);
                clear z1 L3 x v xq;
            end
            clear FN L x1 y1;
        end
        X(1)=[];Y(1)=[];Z(1)=[]; % Deleting first element because horzcat started after 1 initial zero element
        
        %save point_cloud_P3_m3l X Y Z
        Xn{n} = X; Yn{n} = Y; Zn{n} = Z;
        n = n+1;
        clear X Y Z; 

    end
end % SN
clear n muscle MN SN;









