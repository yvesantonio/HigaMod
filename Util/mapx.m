function [newX] = mapx(x,y,map)

    n = length(y);
    m = 
    newX = zeros(1,n);
    
    for ii = 1:n
        mapped = map(x,y(ii));
        newX(ii) = mapped(1);
    end
end