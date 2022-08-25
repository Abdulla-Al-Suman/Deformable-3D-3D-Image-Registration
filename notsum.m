function A = notsum(A,dim)

    %works like a normal sum, but adds in all dimensions except those
    %specified
    y = 1:ndims(A);
    y(dim) = [];
    for i=1:numel(y)
        A = sum(A,y(i));
    end
end