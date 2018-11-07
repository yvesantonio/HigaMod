function [f]=sym2fxy(f)
if(isempty(symvar(f)))
    f=double(f);
    f=@(x,y) f + 0*x+ 0*y;
else
    if(symvar(f)=='x')
        f=matlabFunction(f);
        f=@(x,y) f(x) + 0*y;
    else
        if (symvar(f)=='y')
            f=matlabFunction(f);
            f=@(x,y) f(y) + 0*x;
        else
            f=matlabFunction(f);
        end
    end
end
end