close all;
clear;
tic;

%% Update
% 14 Feb update the damping m.a = F - c/m.v
% 14 Feb test the gravity loading 1D

% Unit
% kN - seconds - m
cfilter=0;

% fid=fopen('elasticSquare1.txt');

%% Material porperties
E                       = 10000     ;                  % Young modulus of solid
psp                     = 4.0           ;                  % solid density
nu                      = 0.41           ;                  % Poison ratio
g                       = 0.0            ;                  % gravity acceleration

D = [E/(1-nu.^2) E*nu/(1-nu.^2) 0; E*nu/(1-nu.^2) E/(1-nu.^2) 0; 0 0 E/(1-nu)];
Lambda                  = E*nu/(1+nu)/(1-2*nu);
Mu                      = E/2/(1+nu);

%% Time
ftime                   = 10.0;
dt                      = 0.001;
ndt                     = round(ftime/dt) +1;
t                       = 0;

% time integration
p_b                     = 1.0;
a_m                     = (2*p_b-1)/(1+p_b);
beta                    = (5-3*p_b)/(1+p_b).^2/(2-p_b);
gamma                   = 3/2-a_m;

%% Grid generation
NN(1,1)                 = 33;                               % number of nodes in X direction
NN(1,2)                 = 13;                               % number of nodes in Y direction
cellCount               = (NN(1)-1)*(NN(2)-1);              % number of elements
nodeCount               = NN(1)*NN(2);                      % number of nodes
LOC                     = zeros(NN(1)*NN(2),2);             % zero matrix of all nodal coordinate
le(1,1)                 = 0.5;                              % size of element in X direction
le(1,2)                 = le(1,1);                          % size of element in Y direction

LOCX                    = (0:NN(1)-1)'*le(1);               % Location of all nodes in X direction
LOCY                    = (0:NN(2)-1)'*le(2);               % Location of all nodes in Y direction

LOCCX                   = (0:(NN(1)-1)-1)'*le(1)+le(1)/2;   % Location of cells in X direction
LOCCY                   = (0:(NN(2)-1)-1)'*le(2)+le(2)/2;   % Location of cells in X direction

for i=1:NN(2)
    LOC((1+NN(1)*(i-1)):(NN(1)*(i-1)+NN(1)),1) = LOCX;      % generate the X node position in LOC
end

for i=1:NN(2)
    LOC((NN(1)*(i-1))+1:NN(1)*i,2) = LOCY(i);               % generate the Y nodeposition in LOC
end

for i=1:NN(2)-1
    LOCC((1+(NN(1)-1)*(i-1)):((NN(1)-1)*(i-1)+(NN(1)-1)),1) = LOCCX;        % generate the X element position in LOCC
end

for i=1:NN(2)-1
    LOCC(((NN(1)-1)*(i-1))+1:(NN(1)-1)*i,2) = LOCCY(i);                     % generate the Y element position in LOCC
end

%% Boundary condition
nfbcx = 0;
nfbcy = 0;
fbcx =[];
fbcy =[];

for n=1:nodeCount
    if LOC(n,1)<=0.5 || LOC(n,1)>=15.5
        nfbcx = nfbcx + 1;
        fbcx = [fbcx n];
    end
end

for n=1:nodeCount
    if LOC(n,2)<=0.5 || LOC(n,2)>=5.5
        nfbcy = nfbcy + 1;
        fbcy = [fbcy n];
    end
end

%% Particle generation
particle_per_cell       = 4;
spCount                 = 100*particle_per_cell;
lp(1,1)                 = le(1,1)/2;                                 % size of particle in X direction
lp(1,2)                 = lp(1,1);                              % size of particle in Y direction
x_sp                    = zeros(spCount,2);
d_sp                    = zeros(spCount,2);

sp=1;
while sp<spCount+0.0001
    for i=1:20
        for j=1:20
            x_sp(sp,1:2)= [1*le(1,1)+0.5*lp(1,1)+(j-1)*lp(1,1) 1*le(1,2)+0.5*lp(1,2)+(i-1)*lp(1,2)];
            sp=sp+1;
        end
    end  
end

sp=1;
while sp<spCount+1
    
        x_center = x_sp(sp,1)-2.5;
        y_center = x_sp(sp,2)-2.5;
        Radius = 1.5;
        
     if (x_center^2+y_center^2)>Radius^2 
    x_sp(sp,:)=[];
    spCount=spCount-1;
    sp=sp-1;
     end
    sp=sp+1;
end
% 
% figure
% plot(x_sp(:,1),x_sp(:,2),'o');
% title('initial condition');
% grid on
% axis([0,16,0,6]);
% set(gca,'xtick',(0:le(1,1):LOCX(NN(1,1))));
% set(gca,'ytick',(0:le(1,2):LOCY(NN(1,2))));
% % set(gca,'XTickLabel',[]);
% % set(gca,'YTickLabel',[]);

%% Initial condition
% Particle
dparticle               = lp(1)*lp(2);                          % area of particle domain
x_spo                   = x_sp;                                 % initial position
b_sp                    = [0 0];                                % body force
s_sp                    = zeros(spCount,3);                     % Stress tensor
ds_sp                   = zeros(spCount,3);                     % Stress increment
v_ssp                   = zeros(spCount,2);                     % velocty
e_sp                    = zeros(spCount,3);                     % Strain tensor
de_sp                   = zeros(spCount,3);                     % Strain increment
ptraction_sp            = zeros(spCount,2);
alpha_sp                = zeros(spCount,1);
r1_sp                   = zeros(spCount,2);
r2_sp                   = zeros(spCount,2);

  mspoints               = cell(cellCount,1);
  F_sp                   = cell(spCount,1);

for sp = 1:spCount
    r1_sp(sp,:) = [lp(1,1)/2 0];
    r2_sp(sp,:) = [0 lp(1,2)/2];
    
       F_sp{sp} = [1 0; 0 1];

end

V_sp                    = zeros(spCount,1);
for sp=1:spCount
V_sp(sp)                = 4*abs(r1_sp(sp,1)*r2_sp(sp,2)-r1_sp(sp,2)*r2_sp(sp,1)); 
end
V_spo                   = V_sp;
m_sp                    = psp * V_sp;                           % mass

for sp=1:spCount
    v_ssp(sp,:) = [5 0];
end

% Node
nmass_si                = zeros(nodeCount,1);
nmomentum_si            = zeros(nodeCount,2);
nvolume_si              = zeros(nodeCount,1);
niforce_si              = zeros(nodeCount,2);
neforce_si              = zeros(nodeCount,2);
nforce_si               = zeros(nodeCount,2);
nvelo_si                = zeros(nodeCount,2);
traction_si             = zeros(nodeCount,2);
nacc_si                 = zeros(nodeCount,2);
nacc_old_si             = zeros(nodeCount,2);
nacc_middle_si          = zeros(nodeCount,2);
nstress_si              = zeros(nodeCount*2,2);
A_si                    = zeros(nodeCount*2,2);

% start the algorithm
% video
timestep = 100;    % number of frame to save
r=timestep/20;      % number of frame per second video ~200s

writerObj2           = VideoWriter('ddmp.avi');
writerObj2.FrameRate = r;    % number of frame per second
open(writerObj2);
    for tt = 1:timestep
         ft              = ftime/timestep*tt;
         
%   ft=ftime;

 while t<ft+0.0000000001      
     t
     
%% Reset value of node
nmass_si(:)                = 0;
nmomentum_si(:)            = 0;  
niforce_si(:)              = 0;
neforce_si(:)              = 0;
nvelo_si(:)                = 0;
nvolume_si(:)              = 0;
nacc_middle_si(:)          = 0;
traction_si(:)             = 0;
nvolume_si(:)              = 0;
nstress_si(:)              = 0;
A_si(:)                    = 0;

%% Store particles into cell

 spElems = zeros(spCount,1);                            % vectors store no of elements for each particle
 CONNECT = zeros(spCount,4);                            % node 1=leftdown 2=righdown 3=rightup 4= leftup
 N = zeros(spCount,4);                                  % Shape function
 dN = zeros(spCount,8);                                 % Gradient of shape function
 
 for sp = 1:spCount
 spElems(sp) = ceil(x_sp(sp,1)/le(1))+(NN(1)-1)*(fix(x_sp(sp,2)/le(2)));   % compute vector store index elements                           
 
 CONNECT(sp,1) = spElems(sp) + floor(spElems(sp)/(NN(1)-1));
 CONNECT(sp,2) = CONNECT(sp,1)+1; 
 CONNECT(sp,3) = CONNECT(sp,2)+NN(1); 
 CONNECT(sp,4) = CONNECT(sp,1)+NN(1);

 [N(sp,1),dN(sp,1),dN(sp,5)]=linearshape(x_sp(sp,1:2),LOC(CONNECT(sp,1),:),le(1,1),le(1,2));
 [N(sp,2),dN(sp,2),dN(sp,6)]=linearshape(x_sp(sp,1:2),LOC(CONNECT(sp,2),:),le(1,1),le(1,2));
 [N(sp,3),dN(sp,3),dN(sp,7)]=linearshape(x_sp(sp,1:2),LOC(CONNECT(sp,3),:),le(1,1),le(1,2));
 [N(sp,4),dN(sp,4),dN(sp,8)]=linearshape(x_sp(sp,1:2),LOC(CONNECT(sp,4),:),le(1,1),le(1,2));
 
 % Build matrix of gradient of shape function of 4 nodes of the cell which
 % contain the particles
 dG((sp-1)*4+(1:4),:) = [dN(sp,1) dN(sp,5);dN(sp,2) dN(sp,6);dN(sp,3) dN(sp,7);dN(sp,4) dN(sp,8)];
 
 % Build alpha
     alpha_sp(sp)=0.5*(256*N(sp,1)*N(sp,2)*N(sp,3)*N(sp,4)).^0.25;

 end
 
 for c =1:cellCount
     id_sp = find(spElems==c);
     mspoints{c}=id_sp;
 end
 
 %% Mapping from particle to nodes
 for sp=1:spCount
 % Build stress tensor
 SSP = [s_sp(sp,1) s_sp(sp,3);s_sp(sp,3) s_sp(sp,2)];
 
 for j=1:4
     npid                           = CONNECT(sp,j);    
     
     if N(sp,j)==0
         continue
     end
     
 % Mass
 nmass_si(npid)            = nmass_si(npid) + m_sp(sp)*N(sp,j);
 
 % Momentum
 nmomentum_si(npid,:)      = nmomentum_si(npid,:) + m_sp(sp)*v_ssp(sp,:)*N(sp,j);
 
 % Internal forces
niforce_si(npid,:)         = niforce_si(npid,:) - (alpha_sp(sp)*V_sp(sp)*SSP*dG((sp-1)*4+j,:)')';

 % External forces
neforce_si(npid,:)         = neforce_si(npid,:) + b_sp*m_sp(sp)*N(sp,j);

% Traction
 traction_si(npid,:)       = traction_si(npid,:) + V_sp(sp)*ptraction_sp(sp,:)*N(sp,j)/le(1,1)/le(1,2);
 
% Volumne
 nvolume_si(npid)          = nvolume_si(npid) + V_sp(sp)*N(sp,j);
 
% Nodal stress
 nstress_si(npid*2-1:npid*2,:)  = nstress_si(npid*2-1:npid*2,:) + (1-alpha_sp(sp))*V_sp(sp)*N(sp,j)*SSP;
end 
 end

 for n=1:nodeCount
     if nvolume_si(n)==0
         continue
     end
     
     nstress_si(n*2-1:n*2,:) = nstress_si(n*2-1:n*2,:)/nvolume_si(n);
 end
 
%% Compute internal force
for c=1:cellCount
    mpts = mspoints{c};
    if isempty(mpts)==1
        continue
    end
    node1 = c + floor(c/(NN(1)-1));
    node2 = node1 + 1;
    node3 = node2 + NN(1);
    node4 = node3 -1;
    niforce_si(node1,:) = niforce_si(node1,:) - [-le(1,2)/6 -le(1,1)/6]     * nstress_si(2*node1-1:2*node1,:);%/nvolume_si(node1);
    niforce_si(node1,:) = niforce_si(node1,:) - [-le(1,2)/6 -le(1,1)/12]    * nstress_si(2*node2-1:2*node2,:);%/nvolume_si(node1);
    niforce_si(node1,:) = niforce_si(node1,:) - [-le(1,2)/12 -le(1,1)/12]   * nstress_si(2*node3-1:2*node3,:);%/nvolume_si(node1);
    niforce_si(node1,:) = niforce_si(node1,:) - [-le(1,2)/12 -le(1,1)/6]    * nstress_si(2*node4-1:2*node4,:);%/nvolume_si(node1);
    
    niforce_si(node2,:) = niforce_si(node2,:) - [le(1,2)/6 -le(1,1)/12]     * nstress_si(2*node1-1:2*node1,:);%/nvolume_si(node2);
    niforce_si(node2,:) = niforce_si(node2,:) - [le(1,2)/6 -le(1,1)/6]      * nstress_si(2*node2-1:2*node2,:);%/nvolume_si(node2);
    niforce_si(node2,:) = niforce_si(node2,:) - [le(1,2)/12 -le(1,1)/6]     * nstress_si(2*node3-1:2*node3,:);%/nvolume_si(node2);
    niforce_si(node2,:) = niforce_si(node2,:) - [le(1,2)/12 -le(1,1)/12]    * nstress_si(2*node4-1:2*node4,:);%/nvolume_si(node2);
    
    niforce_si(node3,:) = niforce_si(node3,:) - [le(1,2)/12 le(1,1)/12]     * nstress_si(2*node1-1:2*node1,:);%/nvolume_si(node3);
    niforce_si(node3,:) = niforce_si(node3,:) - [le(1,2)/12 le(1,1)/6]      * nstress_si(2*node2-1:2*node2,:);%/nvolume_si(node3);
    niforce_si(node3,:) = niforce_si(node3,:) - [le(1,2)/6 le(1,1)/6]       * nstress_si(2*node3-1:2*node3,:);%/nvolume_si(node3);
    niforce_si(node3,:) = niforce_si(node3,:) - [le(1,2)/6 le(1,1)/12]      * nstress_si(2*node4-1:2*node4,:);%/nvolume_si(node3);
    
    niforce_si(node4,:) = niforce_si(node4,:) - [-le(1,2)/12 le(1,1)/6]     * nstress_si(2*node1-1:2*node1,:);%/nvolume_si(node4);
    niforce_si(node4,:) = niforce_si(node4,:) - [-le(1,2)/12 le(1,1)/12]    * nstress_si(2*node2-1:2*node2,:);%/nvolume_si(node4);
    niforce_si(node4,:) = niforce_si(node4,:) - [-le(1,2)/6 le(1,1)/12]     * nstress_si(2*node3-1:2*node3,:);%/nvolume_si(node4);
    niforce_si(node4,:) = niforce_si(node4,:) - [-le(1,2)/6 le(1,1)/6]      * nstress_si(2*node4-1:2*node4,:);%/nvolume_si(node4);
end


%% Update momentum
% Update force
 nforce_si      = niforce_si + neforce_si + traction_si;
 nmomentum_si   = nmomentum_si + nforce_si*dt;
 
 for n=1:nodeCount
     if nmass_si(n)==0
         continue
     end
 nacc_middle_si(n,:) = nforce_si(n,:)/nmass_si(n);
 end
 
 nacc_si           = (nacc_middle_si - nacc_old_si * a_m)/(1 - a_m);
 
  % Boundary condition
% for i=1:nfbc
% nforce_si(fbc(i),:)     = 0;
% nmomentum_si(fbc(i),:)  = 0;
% end

for i=1:nfbcx
nforce_si(fbcx(i),1)     = 0;
nmomentum_si(fbcx(i),1)  = 0;
nvelo_si(fbcx(i),1)  = 0;
end

for i=1:nfbcy
nforce_si(fbcy(i),2)     = 0;
nmomentum_si(fbcy(i),2)  = 0;
nvelo_si(fbcy(i),2)  = 0;
end
 
 %% Update solid particle velocity and position
 % Update position and velocity
  for sp = 1:spCount
      if cfilter == 1
     for j = 1:4
              if nmass_si(CONNECT(sp,j))==0
                continue
              end
     
         v_ssp(sp,:)                      = v_ssp(sp,:) + dt * N(sp,j) * ((1-gamma)*nacc_old_si(CONNECT(sp,j),:)+gamma*nacc_si(CONNECT(sp,j),:));
         x_sp(sp,:)                       = x_sp(sp,:) + N(sp,j)*dt*(nmomentum_si(CONNECT(sp,j),:)/ nmass_si(CONNECT(sp,j)) + dt*((0.5-beta)*nacc_old_si(CONNECT(sp,j),:)+beta*nacc_si(CONNECT(sp,j),:)));
         d_sp(sp,:)                       = d_sp(sp,:) + N(sp,j)*dt*(nmomentum_si(CONNECT(sp,j),:)/ nmass_si(CONNECT(sp,j)) + dt*((0.5-beta)*nacc_old_si(CONNECT(sp,j),:)+beta*nacc_si(CONNECT(sp,j),:)));  
     end
else
     for j = 1:4
              if nmass_si(CONNECT(sp,j))==0
                continue
              end         
         v_ssp(sp,:)                      = v_ssp(sp,:) + dt * N(sp,j) * (nforce_si(CONNECT(sp,j),:)/nmass_si(CONNECT(sp,j)));
         x_sp(sp,:)                       = x_sp(sp,:) + nmomentum_si(CONNECT(sp,j),:)*N(sp,j)*dt/ nmass_si(CONNECT(sp,j));
         d_sp(sp,:)                       = d_sp(sp,:) + nmomentum_si(CONNECT(sp,j),:)*N(sp,j)*dt/ nmass_si(CONNECT(sp,j));
     end   
      end
  end
  
 %% Mapping velocity back to node
 nmomentum_si(:) = 0;
 % Momentum
 for sp = 1:spCount
     for j = 1:4
         if N(sp,j)==0
             continue
         end
         
         % Momentum
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
 
% Boundary condition
% for i=1:nfbc
% nvelo_si(fbc(i),:)     = 0;
% end
  
for i=1:nfbcx
nvelo_si(fbcx(i),1)     = 0;
end

for i=1:nfbcy
nvelo_si(fbcy(i),2)     = 0;
end
  

%% Compute A
for c=1:cellCount
    mpts = mspoints{c};
    if isempty(mpts)==1
        continue
    end

    node1 = c + floor(c/(NN(1)-1));
    node2 = node1 + 1;
    node3 = node2 + NN(1);
    node4 = node3 - 1;
    
    A_si(2*node1-1:2*node1,:) = A_si(2*node1-1:2*node1,:) + nvelo_si(node1,:)'*[-le(1,2)/6 -le(1,1)/6];%/nvolume_si(node1);
    A_si(2*node1-1:2*node1,:) = A_si(2*node1-1:2*node1,:) + nvelo_si(node2,:)'*[le(1,2)/6 -le(1,1)/12];%/nvolume_si(node1);
    A_si(2*node1-1:2*node1,:) = A_si(2*node1-1:2*node1,:) + nvelo_si(node3,:)'*[le(1,2)/12 le(1,1)/12];%/nvolume_si(node1);
    A_si(2*node1-1:2*node1,:) = A_si(2*node1-1:2*node1,:) + nvelo_si(node4,:)'*[-le(1,2)/12 le(1,1)/6];%/nvolume_si(node1);
    
    A_si(2*node2-1:2*node2,:) = A_si(2*node2-1:2*node2,:) + nvelo_si(node1,:)'*[-le(1,2)/6 -le(1,1)/12];%/nvolume_si(node2);
    A_si(2*node2-1:2*node2,:) = A_si(2*node2-1:2*node2,:) + nvelo_si(node2,:)'*[le(1,2)/6 -le(1,1)/6];%/nvolume_si(node2);
    A_si(2*node2-1:2*node2,:) = A_si(2*node2-1:2*node2,:) + nvelo_si(node3,:)'*[le(1,2)/12 le(1,1)/6];%/nvolume_si(node2);
    A_si(2*node2-1:2*node2,:) = A_si(2*node2-1:2*node2,:) + nvelo_si(node4,:)'*[-le(1,2)/12 le(1,1)/12];%/nvolume_si(node2);
    
    A_si(2*node3-1:2*node3,:) = A_si(2*node3-1:2*node3,:) + nvelo_si(node1,:)'*[-le(1,2)/12 -le(1,1)/12];%/nvolume_si(node3);
    A_si(2*node3-1:2*node3,:) = A_si(2*node3-1:2*node3,:) + nvelo_si(node2,:)'*[le(1,2)/12 -le(1,1)/6];%/nvolume_si(node3);
    A_si(2*node3-1:2*node3,:) = A_si(2*node3-1:2*node3,:) + nvelo_si(node3,:)'*[le(1,2)/6 le(1,1)/6];%/nvolume_si(node3);
    A_si(2*node3-1:2*node3,:) = A_si(2*node3-1:2*node3,:) + nvelo_si(node4,:)'*[-le(1,2)/6 le(1,1)/12];%/nvolume_si(node3);
    
    A_si(2*node4-1:2*node4,:) = A_si(2*node4-1:2*node4,:) + nvelo_si(node1,:)'*[-le(1,2)/12 -le(1,1)/6];%/nvolume_si(node4);
    A_si(2*node4-1:2*node4,:) = A_si(2*node4-1:2*node4,:) + nvelo_si(node2,:)'*[le(1,2)/12 -le(1,1)/12];%/nvolume_si(node4);
    A_si(2*node4-1:2*node4,:) = A_si(2*node4-1:2*node4,:) + nvelo_si(node3,:)'*[le(1,2)/6 le(1,1)/12];%/nvolume_si(node4);
    A_si(2*node4-1:2*node4,:) = A_si(2*node4-1:2*node4,:) + nvelo_si(node4,:)'*[-le(1,2)/6 le(1,1)/6];%/nvolume_si(node4);
end

for n=1:nodeCount
     if nvolume_si(n)==0
         A_si(2*n-1:2*n,:)=0;
         continue
     end
     
     A_si(2*n-1:2*n,:)= A_si(2*n-1:2*n,:)/nvolume_si(n);
     
end

%% Update effective stress
for c=1:cellCount
    mpts = mspoints{c};
    dEFSP  = zeros(2,2*length(mpts));
    dEFSP2 = zeros(2,2*length(mpts));
    dS = zeros(2*length(mpts),4);
    
    for sp = 1:length(mpts)
        spid = mpts(sp);
        L_sp = zeros(2,2);
        L_sp2 = zeros(2,2);
        for j=1:4
         if N(sp,j)==0
             continue
         end
            L_sp  = L_sp  + (nvelo_si(CONNECT(spid,j),:)'*dG((spid-1)*4+j,:));
            L_sp2 = L_sp2 + (N(spid,j)*A_si(CONNECT(spid,j)*2-1:CONNECT(spid,j)*2,:));
        end
        
%         Lp_sp(:,sp*2-1:sp*2)    = L_sp;
%         Lp_sp2(:,sp*2-1:sp*2)   = L_sp2;
        dESP                    = (L_sp + L_sp')/2*dt; 
        dEFSP(:,sp*2-1:sp*2)    = dESP;        
        dESP2                   = (L_sp2 + L_sp2')/2*dt; 
        dEFSP2(:,sp*2-1:sp*2)   = dESP2;
        
        L_sp_final = L_sp*alpha_sp(sp) + L_sp2*(1-alpha_sp(sp));
        
        F_sp{spid} = (eye(2,2)+L_sp_final*dt)*F_sp{spid};      
                
        J = det(F_sp{spid});
        V_sp(spid)=V_spo(spid)*J;
         
        SSP = Lambda*log(J)/J*eye(2,2) + Mu/J*(F_sp{spid}*F_sp{spid}'-eye(2,2));
        s_sp(spid,1) = SSP(1,1);
        s_sp(spid,2) = SSP(2,2);
        s_sp(spid,3) = SSP(1,2);
        
        if cfilter==1
        dS(2*sp-1:2*sp,:) = dG((spid-1)*4+1:(spid-1)*4+4,:)';
        end
%     end
    
%     if cfilter==1
%         
%     % QR method
%     if isempty(mpts) ==0
%         [Q1,R]=qr(dS,0);
%         
%         r = rank(R);
% %         if r~=3
% %             break
% %         end
%         
%         dEFSP = (dEFSP*Q1(:,1:r))*Q1(:,1:r)';
%     end    
%     end
    
%     for sp = 1:length(mpts)
%         spid = mpts(sp);
%         de_sp(spid,1) = dEFSP(1,sp*2-1)*alpha_sp(sp)+dEFSP2(1,sp*2-1)*(1-alpha_sp(sp));
%         de_sp(spid,2) = dEFSP(2,sp*2)*alpha_sp(sp)+dEFSP2(2,sp*2)*(1-alpha_sp(sp));
%         de_sp(spid,3) = dEFSP(2,sp*2-1)*alpha_sp(sp)+dEFSP2(2,sp*2-1)*(1-alpha_sp(sp));
%         
%         ds_sp(spid,:) = D * de_sp(spid,:)';
%         s_sp(spid,:)  = s_sp(spid,:) + ds_sp(spid,:);
    end
end
 
 % Update time and step 
 t = t+dt;
 
 % Store old nodal acceleration
 nacc_old_si = nacc_si;
 
 end
 
displacement = zeros(spCount,1);
velocity = zeros(spCount,1);
strain = zeros(spCount,1);
for sp=1:spCount
    displacement(sp) = sqrt(d_sp(sp,1)^2+d_sp(sp,2)^2);
    velocity(sp) = sqrt(v_ssp(sp,1)^2+v_ssp(sp,2)^2);
    strain(sp) = sqrt(e_sp(sp,1)^2+e_sp(sp,2)^2);
end

StressProfile1=figure;
set(StressProfile1, 'visible','off');
sz = 40;
color = strain;
% color = s_sp(:,3);
scatter(x_sp(:,1),x_sp(:,2),sz,color,'filled');
grid on
axis([0,LOCX(NN(1,1)),0,LOCY(NN(1,2))]);
set(gca,'xtick',(0:0.5:LOCX(NN(1,1))));
set(gca,'ytick',(0:0.5:LOCY(NN(1,2))));
colormap(jet(256))

% set(h, 'ytick', [0:0.2:1.2]);
caxis([0 0.25]);

    frame2 = getframe(StressProfile1);
    writeVideo(writerObj2,frame2);
 
    end

    close(writerObj2);