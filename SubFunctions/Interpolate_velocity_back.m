function [nvelo_si] = Interpolate_velocity_back(spCount,CONNECT,nmomentum_si,m_sp,v_ssp,N,nmass_si)

 nmomentum_si(:) = 0;
 
 % Momentum
 for sp = 1:spCount
     for j = 1:4
         nmomentum_si(CONNECT(sp,j),:)  = nmomentum_si(CONNECT(sp,j),:) + m_sp(sp)*v_ssp(sp,:)*N(sp,j);
     end
 end
 
 % Velocity
  for sp = 1:spCount
 for j = 1:4
      if nmass_si(CONNECT(sp,j),:) ==0
         continue
      end
nvelo_si(CONNECT(sp,j),:)               = nmomentum_si(CONNECT(sp,j),:)/nmass_si(CONNECT(sp,j)); 
 end
  end