function [spCount,x_sp,lp]=Two_Particle_Ball(le,x_o,y_o,x_1,y_1,r)

particle_per_cell       = 9;
spCount                 = 22*22*particle_per_cell;
lp(1)                   = le(1)/sqrt(particle_per_cell);      % size of particle in X direction
lp(2)                   = lp(1);                              % size of particle in Y direction
x_sp                    = zeros(spCount,2);

sp=1;
while sp<spCount+0.0001
    for i=1:22*sqrt(particle_per_cell)
        for j=1:22*sqrt(particle_per_cell)
            x_sp(sp,1:2)= [lp(1,1)/2+(j-1)*lp(1,1) lp(1,2)/2+(i-1)*lp(1,2)];
            sp=sp+1;
        end
    end  
end

sp=1;
while sp<spCount+1
        x_center = x_sp(sp,1)-x_o;
        y_center = x_sp(sp,2)-y_o;
        Radius = r;
        x_center1 = x_sp(sp,1)-x_1;
        y_center1 = x_sp(sp,2)-y_1;
        
     if (x_center^2+y_center^2)>Radius^2 && (x_center1^2+y_center1^2)>Radius^2
    x_sp(sp,:)=[];
    spCount=spCount-1;
    sp=sp-1;
     end
    sp=sp+1;
end