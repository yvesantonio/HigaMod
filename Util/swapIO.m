function [mapSwap,JacSwap,HesSwap] = swapIO(geometry)

    aux1 = @(x,y,z) geometry.map([z;y;x]);
    aux2 = @(x,y,z) geometry.map_der([z;y;x]);
    aux3 = @(x,y,z) geometry.map_der2([z;y;x]);
    
    mapSwap(1) = @(x,y,z) aux1(x,y,z);
    mapSwap(2) = aux1(3);
    mapSwap(3) = aux1(1);
    
    JacSwap(1) = aux2(2);
    JacSwap(2) = aux2(3);
    JacSwap(3) = aux2(1);
    
    HesSwap(1) = aux3(2);
    HesSwap(2) = aux3(3);
    HesSwap(3) = aux3(1);

end