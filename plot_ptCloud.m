function plot_ptCloud(X,Y,Z,Xn,Yn,Zn)

for i =1:6
    plot3(X{i},Y{i},Z{i},'.r')
    hold on; plot3(Xn{i},Yn{i},Zn{i},'.k')
    hold on
end
axis equal
hold off
