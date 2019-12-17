function [spCount,x_sp,lp]=Particle_Ball(le)

particle_per_cell       = 9;
spCount                 = 100*particle_per_cell;
lp(1)                   = le(1)/sqrt(particle_per_cell);                            % size of particle in X direction
lp(2)                   = lp(1);                              % size of particle in Y direction
x_sp                    = zeros(spCount,2);

sp=1;
while sp<spCount+0.0001
    for i=1:10*sqrt(particle_per_cell)
        for j=1:10*sqrt(particle_per_cell)
            x_sp(sp,1:2)= [1*le(1)+0.5*lp(1,1)+(j-1)*lp(1,1) 1*le(2)+0.5*lp(1,2)+(i-1)*lp(1,2)];
            sp=sp+1;
        end
    end  
end

sp=1;
while sp<spCount+1
        x_center = x_sp(sp,1)-2.5;
        y_center = x_sp(sp,2)-3;
        Radius = 1.5;
        
     if (x_center^2+y_center^2)>Radius^2 
    x_sp(sp,:)=[];
    spCount=spCount-1;
    sp=sp-1;
     end
    sp=sp+1;
end