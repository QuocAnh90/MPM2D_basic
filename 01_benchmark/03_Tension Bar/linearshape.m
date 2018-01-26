function [N,dN1,dN2]=linearshape(xp,xn,Lx,Ly)


if abs(xp(1)-xn(1))<Lx 
    
    Nx = 1-abs(xp(1)-xn(1))/Lx;    
    dNx = -sign(xp(1)-xn(1))/Lx;
else
    Nx = 0;
    dNx = 0;
end

if abs(xp(2)-xn(2))<Ly 
    Ny = 1-abs(xp(2)-xn(2))/Ly;     
    dNy = -sign(xp(2)-xn(2))/Ly;
else
    Ny = 0;
    dNy = 0;
end


    N = Nx*Ny;
    dN1 = dNx*Ny;
    dN2 = Nx*dNy;
   
end