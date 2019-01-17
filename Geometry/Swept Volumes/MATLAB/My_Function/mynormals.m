function [Nx Ny Nz]=mynormals(p1,p2,p3)
    u=p2-p1;
    v=p3-p1;
    Nx=u(2)*v(3)-u(3)*v(2);
    Ny=u(1)*v(1)-u(1)*v(3);
    Nz=u(1)*v(2)-u(2)*v(1);
end