function  bx=NeumannInflow(wyq,bx,Pinf,MbRX,fedof,J)

for i=1:size(MbRX,2)
   bx(1+ (i-1)*fedof) = bx(1 +(i-1)*fedof)+Pinf*J(:,1)'*(MbRX(:,i).*wyq);
end
end