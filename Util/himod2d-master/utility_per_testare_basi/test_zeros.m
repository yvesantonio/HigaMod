addpath ../
global mu coeffrobin
mu=@(x,y) 1;
coeffrobin = 12;
muval=mu(0,0);
sigma = coeffrobin;
fplot(@(lambda) 2*muval*lambda + tan(lambda)*(sigma - muval*lambda.^2/sigma),[0,30]);
hold on
fplot(@(x) 0,[0,30])
ylim([-20,20])
