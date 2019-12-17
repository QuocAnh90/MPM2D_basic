function [v_ssp,x_sp,d_sp] = Update_Particle_Position(NODES,dt,CONNECT,N,spCount,nmass_si,nforce_si,nmomentum_si,x_spo,v_ssp,x_sp,d_sp)

for sp = 1:spCount
     for j = 1:NODES(sp)
         npid                           = CONNECT{sp}(j);
              if nmass_si(npid)==0
                continue
              end         
         v_ssp(sp,:)                      = v_ssp(sp,:) + dt * N{sp}(j) * (nforce_si(npid,:)/nmass_si(npid));
         x_sp(sp,:)                       = x_sp(sp,:) + nmomentum_si(npid,:)*N{sp}(j)*dt/ nmass_si(npid);
         d_sp(sp,:)                       = x_sp(sp,:) - x_spo(sp,:);
     end   
end