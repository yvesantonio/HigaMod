function  bx=NeumannOutflow(wyq,bx,Pout,MbRX,fedof,J)

for i=1:size(MbRX,2)
    bx( i*fedof) = bx( i*fedof )+Pout*J(:,2)'*(MbRX(:,i).*wyq);
end
end
