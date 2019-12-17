function [N,dNx,dNy]=cubic_Bsplineshape(xp,xn,Lx,Ly)       

dx = xp(1)-xn(1);

if dx >= -2*Lx && dx <= -Lx
    Nx = 1/6/Lx^3*dx^3+1/Lx^2*dx^2+2*dx/Lx+4/3;
    dNx = 1/2/Lx^3*dx^2+2/Lx^2*dx+2/Lx;
    
elseif dx > -Lx && dx <= 0
    Nx = -1/2/Lx^3*dx^3-1/Lx^2*dx^2+2/3;
    dNx = -3/2/Lx^3*dx^2-2/Lx^2*dx;
    
elseif dx > 0 && dx <= Lx
    Nx = 1/2/Lx^3*dx^3-1/Lx^2*dx^2+2/3;
    dNx = 3/2/Lx^3*dx^2-2/Lx^2*dx;
    
elseif dx > Lx && dx <= 2*Lx
    Nx = -1/6/Lx^3*dx^3+1/Lx^2*dx^2-2*dx/Lx+4/3;
    dNx = -1/2/Lx^3*dx^2+2/Lx^2*dx-2/Lx;
    
else
    Nx = 0;
    dNx = 0;
end

% dy = xp(2)-xn(2);
% 
% if dy >= -2*Ly && dy <= -Ly
%     Ny = 1/6/Ly^3*dy^3+1/Ly^2*dy^2+2*dy/Ly+4/3;
%     dNy = 1/2/Ly^3*dy^2+2/Ly^2*dy+2/Ly;
%     
% elseif dy > -Ly && dy <= 0
%     Ny = -1/2/Ly^3*dy^3-1/Ly^2*dy^2+2/3;
%     dNy = -3/2/Ly^3*dy^2-2/Ly^2*dy;
%     
% elseif dy > 0 && dy <= Ly
%     Ny = 1/2/Ly^3*dy^3-1/Ly^2*dy^2+2/3;
%     dNy = 3/2/Ly^3*dy^2-2/Ly^2*dy;
%     
% elseif dy > Ly && dy <= 2*Ly
%     Ny = -1/6/Ly^3*dy^3+1/Ly^2*dy^2-2*dy/Ly+4/3;
%     dNy = -1/2/Ly^3*dy^2+2/Ly^2*dy-2/Ly;
%     
% else
%     Ny = 0;
%     dNy = 0;
% end

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
