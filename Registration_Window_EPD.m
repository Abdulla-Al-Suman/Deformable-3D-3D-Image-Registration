function [R_Moving_V_WindowEPD,R_Parameters_WindowEPD] = Registration_Window_EPD(Fixed_Volume,Moving_Volume,d1,d2,d3,LT,UT,x,y,z)
% global maxChamp lev
%global row col
% Creating chamfer distance od fixed and edgelist of moving volumes
EF = zeros(d1,d2,d3);
for i=1:d3
EF(:,:,i) = edge(Fixed_Volume(:,:,i),'canny',[LT UT],1.5); % Calculation of Canny Edge detection of moving image volume
end

DF = ChamferDis3dinfv2(uint32(EF),inf);
clear EF i;
[Edge_list_Moving,Number_EdgePoints_Moving] = Calculate_EdgeList(Moving_Volume,d1,d2,d3,LT,UT);

%maxChamp = inf;

%D = ChamferDis2D(uint32(Er),double(maxChamp),row,col);

% convert binary edge map into list of pixel positions
%edge_list = Bin2List(double(Ei));


% rotation around the point (Cx,Cy)
% Cx = median(edge_list(:,2));
% Cy = median(edge_list(:,1));
% Cx = 512;
% Cy = 0;

% create 4 x 1 array of non-zero pixel positions
%edge_list
%[n,~] = size(edge_list);
v = [Edge_list_Moving ones(Number_EdgePoints_Moving,1)];
clear Edge_list_Moving Number_EdgePoints_Moving;

% initialize EPD similarity measure
Smin = inf;

% w = 6; % search region -w to w

% conduct full search for in-plane position which produces minimum EPD
for tx = -3:1:3
    for ty = -4:1:4 %-w-20:4:w+20
%       M = makeM2d(-tx,-ty,0,0,-Rz,0,0,Cy,Cx);
        M = makeM3d(tx,ty,0,1,1,1);
        
        % calculate new list of edge positions
        vn = v*M';
        new_edge_list=vn(:,1:3);
        new_edge_list=round(new_edge_list);
        new_edge_list(new_edge_list(:, 1)>d2, :)= [];
        new_edge_list(new_edge_list(:, 1)<1, :)= [];
        new_edge_list(new_edge_list(:, 2)>d1, :)= [];
        new_edge_list(new_edge_list(:, 2)<1, :)= [];
        new_edge_list(new_edge_list(:, 3)>d3, :)= [];
        new_edge_list(new_edge_list(:, 3)<1, :)= [];
        % find integer and fractional components of edge positions
        %[xi,xf] = IntFrac(vn);
        
        %xi = uint32(xi);
        
        % calculate the EPD similarity measure at DRR edge positions
        %S = CalcEPD(xi,xf,D);
        S = CalcEPD3D(new_edge_list,DF);
        
        % update minimum EPD value and in-plane parameters if EPD value is
        % less than previous minimum EPD
        if S < Smin
            Smin = S;
            txm = tx;
            tym = ty;
        end
        clear M vn new_edge_list ty;
    end
    clear tx;
end
R_Parameters_WindowEPD=[1.0000 txm 1.0000 tym 1.0000 0.0000];
R_Moving_V_WindowEPD = warp_Window_ST_3D(Moving_Volume,R_Parameters_WindowEPD,x,y,z);
clear Smin txm tym;


% M = makeM2d(-txm,-tym,0,0,-rzm.*pi/180,0,0,Cy,Cx);
% vn = v*M';	
% edge_list = round(vn(:,1:2));
% 
% % delete rows contains less than zero or greater than length of the row
% edge_list(any(edge_list<1,2),:) = [];
% edge_list(any(edge_list>row,2),:) = [];
% 
% % m = -[rzm.*pi/180,txm,tym];
% m = [-rzm.*pi/180 tym txm Cx-256.5 Cy-256.5];

