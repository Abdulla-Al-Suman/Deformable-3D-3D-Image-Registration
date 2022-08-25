function [Edge_list,Number_EdgePoints] = Calculate_EdgeList(Fixed_Volume,d1,d2,d3,LT,UT)


ER = zeros(d1,d2,d3);
for i=1:d3
    ER(:,:,i) = edge(Fixed_Volume(:,:,i),'canny',[LT UT],1.5);
end
[r,c,v] = ind2sub(size(ER),find(ER == 1)); % convert binary edge map into list of pixel positions
Edge_list = [c r v];
clear r c v ER i;
[Number_EdgePoints,~] = size(Edge_list);