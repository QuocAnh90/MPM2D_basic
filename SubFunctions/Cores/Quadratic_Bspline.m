function [N]=Quadratic_Bspline(xp,xn,Lx,Ly)

dx = -xp(1)+xn(1);
rx = abs(dx)/Lx;

if (rx>1.5)
     Nx     = 0.0;
    elseif (rx<=0.5)
     Nx     = 3/4 - rx^2;
    else
     Nx     = 0.5 * (1.5 - rx)^2;
end
    
dy = -xp(2)+xn(2);
ry = abs(dy)/Ly;

if (ry>1.5)
     Ny     = 0.0;
    elseif (ry<=0.5)
     Ny     = 3/4 - ry^2;
    else
     Ny     = 0.5 * (1.5 - ry)^2;
end

    N = Nx*Ny;  
end
