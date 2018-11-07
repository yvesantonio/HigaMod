function f=sym2fxt(f)
    if(isempty(symvar(f)))
        f=double(f);
        f=@(t,x) f + 0*x + 0*t;
    else
        if(symvar(f)=='x')
            f=matlabFunction(f);    
            f=@(t,x) f(x) + 0*t;
        else
            if(symvar(f)=='t')
                f=matlabFunction(f);    
                f=@(t,x) f(t) + 0*x;
            else
                f=matlabFunction(f);
            end
        end
    end
end