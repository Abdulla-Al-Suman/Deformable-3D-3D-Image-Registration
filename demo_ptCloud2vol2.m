clc
clear
close all

load('transfered_points_3F_13M.mat')

load('Patient 13CI muscle data.mat')
[Xm,Ym,Zm] = get_ptCloud(muscle_data);

Vol = ptCloud2vol(Xm,Ym,Zm,3);

for k = 1:128
    imagesc(Vol(:,:,k))
    k
     pause
end