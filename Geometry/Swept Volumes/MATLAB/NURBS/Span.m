function i=Span(n,p,u,U)

if u==U(n)
    i=n;
else
if u==U(p)
        i=p;
    else
        low=p;up=n;mid=fix((low+up)/2); % l lower , u  upper , m mid
        while (u<U(mid) | u>=U(mid+1))
            if u<U(mid)
                up=mid;
            else
                low=mid;
            end
            mid=ceil((low+up)/2);
        end
        i=mid
end
end

end

            