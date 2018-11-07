function  bx=NeumannOutflowREL(wyq,bx,Pout,MbRX,fedof,J)

for i=3:size(MbRX,2)
    %when the domain will be moving check Map.Jac should not be end, but a
    %different one!!!!
    bx( (i-2)*fedof+2) = bx( (i-2)*fedof +2)-Pout*J(:,2)'*(MbRX(:,i).*wyq);
end
end
