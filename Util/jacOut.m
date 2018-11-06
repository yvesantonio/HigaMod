function [evalJac] = jacOut(x,y,Jac)

    n = length(y);
    m = length(x);
    
    for ii = 1:n
        for jj = 1:m
            evalJac{ii,jj} = Jac(x(jj),y(ii));
        end
    end

end