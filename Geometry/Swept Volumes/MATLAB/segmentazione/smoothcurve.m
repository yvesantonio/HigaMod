function St=smoothcurve(St)
if St.cont==1
    St.next1=smoothcurve(St.next1);
    St.next2=smoothcurve(St.next2);
end
[St.data(:,1:3)]=[nrb_approx(St.data(:,1),St.data(:,2),St.data(:,3),max(3,round(numel(St.data(:,1))/10)),2,1.e-5,10,20)]';
% =[nrbeval(N,linspace(0,1,numel(St.data(:,1))))]';
disp('aaaa')