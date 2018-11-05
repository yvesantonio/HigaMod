function [selectedOut] = invJac(ith,jth,evalJac)

    [n,m] = size(evalJac);
    selectedOut = [];
    
    for ii = 1:n
        for jj = 1:m
            evalJavInv = evalJac{ii,jj}\eye(size(evalJac{ii,jj}));
            selectedOut(ii,jj) = evalJavInv(ith,jth);
        end
    end

end