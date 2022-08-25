function M = makeM3d(Tx,Ty,Tz,Sx,Sy,Sz)
    
% Construct translation matrix

Mt = [1 0 0 Tx;...
    0 1 0 Ty;...
    0 0 1 Tz;...
    0 0 0 1];

% Construct zoom matrix

Ms = [Sx 0 0 0;...
    0 Sy 0 0;...
    0 0 Sz 0;...
    0 0 0 1];

% Construct final M matrix

M = Mt*Ms;
% M = Mt*Ma*Mrz*Mb;
clear Mt Ms;