function [N,dN1,dN2]=GIMPshape(xp,xn,Lx,Ly,dparticle)   

% dparticle: dimension of 1 size of particle

lp = dparticle/2;
dx = xp(1) - xn(1);
dy = xp(2) - xn(2);

if dx >= -Lx-lp && dx <= -Lx+lp
    Nx = (Lx+lp+dx).^2/(4*Lx*lp);
    dNx = (Lx+lp+dx)/(2*Lx*lp);
    
elseif dx > -Lx+lp && dx <= -lp
    Nx = 1 + dx/Lx;
    dNx = 1/Lx;
    
elseif dx > -lp && dx <= lp
    Nx = 1 - (dx.^2+lp.^2)/(2*Lx*lp);
    dNx = -dx/Lx/lp;
    
elseif dx > lp && dx <= Lx-lp
    Nx = 1-dx/Lx;
    dNx = -1/Lx;
    
elseif dx > Lx-lp && dx <= Lx+lp
    Nx = (Lx+lp-dx).^2/(4*Lx*lp);
    dNx = -(Lx+lp-dx)/(2*Lx*lp);
    
else
    Nx = 0;
    dNx = 0;
end 

if dy >= -Ly-lp && dy <= -Ly+lp
    Ny = (Ly+lp+dy).^2/(4*Ly*lp);
    dNy = (Ly+lp+dy)/(2*Ly*lp);
    
elseif dy > -Ly+lp && dy <= -lp
    Ny = 1 + dy/Ly;
    dNy = 1/Ly;
    
elseif dy > -lp && dy <= lp
    Ny = 1 - (dy.^2+lp.^2)/(2*Ly*lp);
    dNy = -dy/Ly/lp;
    
elseif dy > lp && dy <= Ly-lp
    Ny = 1-dy/Ly;
    dNy = -1/Ly;
    
elseif dy > Ly-lp && dy <= Ly+lp
    Ny = (Ly+lp-dy).^2/(4*Ly*lp);
    dNy = -(Ly+lp-dy)/(2*Ly*lp);
    
else
    Ny = 0;
    dNy = 0;
end

    N = Nx*Ny;
    dN1 = dNx*Ny;
    dN2 = Nx*dNy;

    end