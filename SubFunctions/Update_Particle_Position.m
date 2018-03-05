function [v_ssp,x_sp,d_sp] = Update_Particle_Position(dt,CONNECT,N,spCount,nmass_si,nforce_si,nmomentum_si,x_spo,v_ssp,x_sp)

for sp = 1:spCount
     for j = 1:4
              if nmass_si(CONNECT(sp,j))==0
                continue
              end         
         v_ssp(sp,:)                      = v_ssp(sp,:) + dt * N(sp,j) * (nforce_si(CONNECT(sp,j),:)/nmass_si(CONNECT(sp,j)));
         x_sp(sp,:)                       = x_sp(sp,:) + nmomentum_si(CONNECT(sp,j),:)*N(sp,j)*dt/ nmass_si(CONNECT(sp,j));
         d_sp(sp,:)                       = x_sp(sp,:) - x_spo(sp,:);
     end
end