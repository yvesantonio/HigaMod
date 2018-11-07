function  bx=NeumannInflowREL(wyq,bx,Pinf,MbRX,fedof,J)

for i=1:size(MbRX,2)
    if(i>2)
    bx(1 +2 + (i-1-2)*fedof) = bx(1 +2 +(i-1-2)*fedof)+Pinf*J(:,1)'*(MbRX(:,i).*wyq);
    else
        
    bx(i) = bx(i)+Pinf*J(:,1)'*(MbRX(:,i).*wyq);
    end
end
end