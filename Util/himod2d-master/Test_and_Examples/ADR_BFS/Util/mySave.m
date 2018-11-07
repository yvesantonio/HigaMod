function mySave(XR,XL,YR,YL)
N=size(XR,1)*size(XR,2)+size(XL,1)*size(XL,2);
file =fopen('./tmp/grid.dat','w');
fprintf(file,'%i\n',N);
for k=1:size(XR,1) %y
    for j=1:size(XR,2) %x
        fprintf(file,'%12.8d     %12.8f\n',[XR(k,j);YR(k,j)]);
    end
end
for k=1:size(XL,1)
    for j=1:size(XL,2)
        fprintf(file,'%12.8d     %12.8f\n',[XL(k,j);YL(k,j)]);
    end
end
fclose(file);