function [evalDetJac] = detJac(x,y,Jac)

    n = length(y);
    m = length(x);
    
    for ii = 1:n
        for jj = 1:m
            evalDetJac(ii,jj)  = abs(det(Jac(x(jj),y(ii))));
        end
    end

end