function [f]=sym2fxyt(f)
if(isempty(symvar(f)))
    f=double(f);
    f=@(t,x,y) f + 0*x+ 0*y + 0*t;
else
    if(length(symvar(f))==1)
        %if(symvar(f)=='x')
        if( isequal(symvar(f),sym('x')))
            f=matlabFunction(f);
            f=@(t,x,y) f(x) + 0*y + 0*t;
            %elseif (symvar(f)=='y')
        elseif( isequal(symvar(f),sym('y')))
            f=matlabFunction(f);
            f=@(t,x,y) f(y) + 0*x + 0*t;
            
            %elseif (symvar(f)=='t')
        elseif( isequal(symvar(f),sym('t')))
            f=matlabFunction(f);
            f=@(t,x,y) f(t) + 0*x + 0*y;
        end
    else
        if(length(symvar(f))==2)
            %if (sum(symvar(f)==[sym('x'),sym('y')])==2)
            if( isequal(symvar(f),[sym('x'),sym('y')]) )
                f=matlabFunction(f);
                f=@(t,x,y) f(x,y) + 0*t;
            else
                %if (sum(symvar(f)==[sym('t'),sym('x')])==2)
                if( isequal(symvar(f),[sym('t'),sym('x')]) )
                    f=matlabFunction(f);
                    f=@(t,x,y) f(t,x)+0*y;
                else
                    %if (sum(symvar(f)==[sym('t'),sym('y')])==2)
                    if( isequal(symvar(f),[sym('t'),sym('y')]) )
                        f=matlabFunction(f);
                        f=@(t,x,y) f(t,y)+0*x;
                    end
                end
            end
        else
            f=matlabFunction(f);
        end
    end
end
end