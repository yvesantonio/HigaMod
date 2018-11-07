function f=sym2fx(f)
    if(isempty(symvar(f)))
        f=@(x) f + 0*x;
    else
        f=matlabFunction(f);
    end
end