function [b]=thereis(v,i)
toll=1e-10;
for k=1:length(v)
    if(abs(i-v(k))<toll)
        b=true;
        return;
    end
end
b=false;
return;