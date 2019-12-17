function [N]=Cubic_Bspline(xp,xn,Lx,Ly)

dx = -xp(1)+xn(1);
rx = abs(dx)/2/Lx;

if (rx>1.0)
     Nx     = 0.0;
    elseif (rx<=0.5)
     Nx     = 2/3 - 4*rx^2 + 4*rx^3;
    else
     Nx     = 4/3 - 4*rx + 4*rx^2 - 4*rx^3/3;
end
    
dy = -xp(2)+xn(2);
ry = abs(dy)/2/Ly;

if (ry>1.0)
     Ny     = 0.0;
    elseif (ry<=0.5)
     Ny     = 2/3 - 4*ry^2 + 4*ry^3;
    else
     Ny     = 4/3 - 4*ry + 4*ry^2 - 4*ry^3/3;
end

    N = Nx*Ny;  
end
