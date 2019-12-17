function [nvelo] = Interpolate_velocity_back(NODES,nodeCount,pCount,CONNECT,m_p,v_p,N,nmass)

%% Interpolation from particle to grid task
% Node variables
 nvelo                = zeros(nodeCount,2);                   % Nodal Velocity
 nmomentum_si            = zeros(nodeCount,2);                   % Nodal Momentum

 % Momentum
 for sp = 1:pCount
     for j = 1:NODES(sp)
         npid                  = CONNECT{sp}(j);
         nmomentum_si(npid,:)  = nmomentum_si(npid,:) + m_p(sp)*v_p(sp,:)*N{sp}(j);
     end
 end
 
 % Velocity
  for sp = 1:pCount
 for j = 1:NODES(sp)
     npid                      = CONNECT{sp}(j);
              if nmass(npid)==0
                continue
              end 
              
nvelo(npid,:)               = nmomentum_si(npid,:)/nmass(npid); 
 end
  end