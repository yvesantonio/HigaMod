function [newOut] = mapOut(x,y,map,out)

    n = length(y);
    m = length(x);
    newOut = zeros(n,m);

    for ii = 1:n
        for jj = 1:m            
            mapped = map(x(jj),y(ii));
            newOut(ii,jj) = mapped(out);
        end
    end
    
end