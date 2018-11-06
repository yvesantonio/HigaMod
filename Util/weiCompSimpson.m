function [w] = weiCompSimpson(n,a,b)

    w = zeros(n+1,1);
    h = (b - a)/n;
    
    w(1)   = h/3;
    w(end) = h/3;
    
    for ii =1:(n/2)
        w(2 * ii) = 4*h/3;
    end
    for ii =1:(n/2 - 1)
        w(2 * ii + 1) = 2*h/3;
    end

end