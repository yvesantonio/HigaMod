function test_ortonormal(f,g,a,b)
h=@(x) f(x).*g(x);
quad(h,a,b)
return 