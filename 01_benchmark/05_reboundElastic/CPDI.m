% close all;
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
ftime                   = 10;
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

figure
plot(x_sp(:,1),x_sp(:,2),'o');
title('initial condition');
grid on
axis([0,16,0,6]);
set(gca,'xtick',(0:le(1,1):LOCX(NN(1,1))));
set(gca,'ytick',(0:le(1,2):LOCY(NN(1,2))));
% set(gca,'XTickLabel',[]);
% set(gca,'YTickLabel',[]);

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

 spElems                = zeros(spCount,1);                     % List of cell index for each particle
 spElems_corner         = zeros(spCount,4);                     % List of cell index for each particle corners
 CONNECT_TEMP           = zeros(spCount,16);                    % List of node index for each particle corners
 NODES                  = zeros(spCount,1);                     % Number of interacting nodes
 CONNECT                = cell(spCount,1);
 N                      = cell(spCount,1);
 dN                     = cell(spCount,1);
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

r10_sp = r1_sp;
r20_sp = r2_sp;

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

%% start the algorithm
% video
timestep = 100;    % number of frame to save
r=timestep/20;      % number of frame per second video ~200s

writerObj2           = VideoWriter('cpdi.avi');
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
% nforce_si(:)               = 0;
niforce_si(:)              = 0;
neforce_si(:)              = 0;
nvelo_si(:)                = 0;
nvolume_si(:)              = 0;
% nacc_si(:)                 = 0;
nacc_middle_si(:)          = 0;
traction_si(:)             = 0;

%% Store particles into cell

 x_corner1 = zeros(spCount,2);
 x_corner2 = zeros(spCount,2);
 x_corner3 = zeros(spCount,2);
 x_corner4 = zeros(spCount,2);
 
 for sp=1:spCount
 x_corner1(sp,:) = x_sp(sp,:) - r1_sp(sp,:) - r2_sp(sp,:);      % Position of corner 1
 x_corner2(sp,:) = x_sp(sp,:) + r1_sp(sp,:) - r2_sp(sp,:);
 x_corner3(sp,:) = x_sp(sp,:) + r1_sp(sp,:) + r2_sp(sp,:);
 x_corner4(sp,:) = x_sp(sp,:) - r1_sp(sp,:) + r2_sp(sp,:);
 end
 
 for sp = 1:spCount
 spElems(sp) = ceil(x_sp(sp,1)/le(1))+(NN(1)-1)*(fix(x_sp(sp,2)/le(2)));   % compute vector store index elements                           

 spElems_corner(sp,1) = ceil(x_corner1(sp,1)/le(1))+(NN(1)-1)*(fix(x_corner1(sp,2)/le(2)));                        
 spElems_corner(sp,2) = ceil(x_corner2(sp,1)/le(1))+(NN(1)-1)*(fix(x_corner2(sp,2)/le(2)));                        
 spElems_corner(sp,3) = ceil(x_corner3(sp,1)/le(1))+(NN(1)-1)*(fix(x_corner3(sp,2)/le(2)));
 spElems_corner(sp,4) = ceil(x_corner4(sp,1)/le(1))+(NN(1)-1)*(fix(x_corner4(sp,2)/le(2)));

 CONNECT_TEMP(sp,1)  = spElems_corner(sp,1) + floor(spElems_corner(sp,1)/(NN(1)-1));
 CONNECT_TEMP(sp,2)  = CONNECT_TEMP(sp,1) + 1;
 CONNECT_TEMP(sp,3)  = CONNECT_TEMP(sp,2) + NN(1);
 CONNECT_TEMP(sp,4)  = CONNECT_TEMP(sp,1) + NN(1);
 
 CONNECT_TEMP(sp,5)  = spElems_corner(sp,2) + floor(spElems_corner(sp,2)/(NN(1)-1));
 CONNECT_TEMP(sp,6)  = CONNECT_TEMP(sp,5) + 1;
 CONNECT_TEMP(sp,7)  = CONNECT_TEMP(sp,6) + NN(1);
 CONNECT_TEMP(sp,8)  = CONNECT_TEMP(sp,5) + NN(1); 

 CONNECT_TEMP(sp,9)  = spElems_corner(sp,3) + floor(spElems_corner(sp,3)/(NN(1)-1));
 CONNECT_TEMP(sp,10) = CONNECT_TEMP(sp,9) + 1;
 CONNECT_TEMP(sp,11) = CONNECT_TEMP(sp,10) + NN(1);
 CONNECT_TEMP(sp,12) = CONNECT_TEMP(sp,9) + NN(1); 

 CONNECT_TEMP(sp,13) = spElems_corner(sp,4) + floor(spElems_corner(sp,4)/(NN(1)-1));
 CONNECT_TEMP(sp,14) = CONNECT_TEMP(sp,13) + 1;
 CONNECT_TEMP(sp,15) = CONNECT_TEMP(sp,14) + NN(1);
 CONNECT_TEMP(sp,16) = CONNECT_TEMP(sp,13) + NN(1);

 CONNECT{sp}=unique(CONNECT_TEMP(sp,:));    % Store nodes interacting with corners
 NODES(sp)=length(CONNECT{sp});         % Store number of interacting nodes
 
 N1 = zeros(1,NODES(sp));   % Shape function of corner 1 for 16 nodes
 N2 = zeros(1,NODES(sp));
 N3 = zeros(1,NODES(sp));
 N4 = zeros(1,NODES(sp));
 N_local = zeros(1,NODES(sp));
 
 for i=1:NODES(sp)
    [N1(i),~,~] = linearshape(x_corner1(sp,:),LOC(CONNECT{sp}(i),:),le(1,1),le(1,2));
    [N2(i),~,~] = linearshape(x_corner2(sp,:),LOC(CONNECT{sp}(i),:),le(1,1),le(1,2));
    [N3(i),~,~] = linearshape(x_corner3(sp,:),LOC(CONNECT{sp}(i),:),le(1,1),le(1,2));
    [N4(i),~,~] = linearshape(x_corner4(sp,:),LOC(CONNECT{sp}(i),:),le(1,1),le(1,2));
    N_local(i)  = 0.25*(N1(i)+N2(i)+N3(i)+N4(i));
 end
 
 N{sp} = N_local;
 
 % Build matrix of gradient of shape function 
 w1 = [r1_sp(sp,2)-r2_sp(sp,2) r2_sp(sp,1)-r1_sp(sp,1)];
 w2 = [r1_sp(sp,2)+r2_sp(sp,2) -r2_sp(sp,1)-r1_sp(sp,1)];
 
 for i=1:NODES(sp)      
dN{sp}(1,i)     = 1/V_sp(sp)*((N1(i)-N3(i))*w1(1) + (N2(i)-N4(i))*w2(1));
dN{sp}(2,i)     = 1/V_sp(sp)*((N1(i)-N3(i))*w1(2) + (N2(i)-N4(i))*w2(2));
 end
 end
 
 for c =1:cellCount
     id_sp = find(spElems==c);
     mspoints{c}=id_sp;
 end
 
 %% Mapping from particle to nodes
 for sp=1:spCount
 % Build stress tensor
 SSP = [s_sp(sp,1) s_sp(sp,3);s_sp(sp,3) s_sp(sp,2)];
 
 for j=1:NODES(sp)
     npid                           = CONNECT{sp}(j);
     
          if N{sp}(j)==0
         continue
          end
     
 % Mass
 nmass_si(npid)            = nmass_si(npid) + m_sp(sp)*N{sp}(j);
 
 % Momentumlp(1,1)lp(1,1)lp(1,1)
 nmomentum_si(npid,:)      = nmomentum_si(npid,:) + m_sp(sp)*v_ssp(sp,:)*N{sp}(j);
 
 % Internal forces
niforce_si(npid,:)         = niforce_si(npid,:) - (V_sp(sp)*SSP*dN{sp}(:,j))';

 % External forces
neforce_si(npid,:)         = neforce_si(npid,:) + b_sp*m_sp(sp)*N{sp}(j);

% Traction
 traction_si(npid,:)       = traction_si(npid,:) + V_sp(sp)*ptraction_sp(sp,:)*N{sp}(j)/le(1,1)/le(1,2);
 
end 
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
     for j = 1:NODES(sp)
         npid                           = CONNECT{sp}(j);
              if nmass_si(npid)==0
                continue
              end
     
         v_ssp(sp,:)                      = v_ssp(sp,:) + dt * N{sp}(j) * ((1-gamma)*nacc_old_si(npid,:)+gamma*nacc_si(npid,:));
         x_sp(sp,:)                       = x_sp(sp,:) + N{sp}(j)*dt*(nmomentum_si(npid,:)/ nmass_si(npid) + dt*((0.5-beta)*nacc_old_si(npid,:)+beta*nacc_si(npid,:)));
         d_sp(sp,:)                       = d_sp(sp,:) + N{sp}(j)*dt*(nmomentum_si(npid,:)/ nmass_si(npid) + dt*((0.5-beta)*nacc_old_si(npid,:)+beta*nacc_si(npid,:)));  
     end
else
     for j = 1:NODES(sp)
         npid                           = CONNECT{sp}(j);
              if nmass_si(npid)==0
                continue
              end         
         v_ssp(sp,:)                      = v_ssp(sp,:) + dt * N{sp}(j) * (nforce_si(npid,:)/nmass_si(npid));
         x_sp(sp,:)                       = x_sp(sp,:) + nmomentum_si(npid,:)*N{sp}(j)*dt/ nmass_si(npid);
         d_sp(sp,:)                       = d_sp(sp,:) + nmomentum_si(npid,:)*N{sp}(j)*dt/ nmass_si(npid);
     end   
end
 end
 
 %% Mapping velocity back to node
 nmomentum_si(:) = 0;
 % Momentum
 for sp = 1:spCount
     for j = 1:NODES(sp)
         npid                  = CONNECT{sp}(j);
         nmomentum_si(npid,:)  = nmomentum_si(npid,:) + m_sp(sp)*v_ssp(sp,:)*N{sp}(j);
     end
 end
 
 % Velocity
  for sp = 1:spCount
 for j = 1:NODES(sp)
     npid                      = CONNECT{sp}(j);
              if nmass_si(npid)==0
                continue
              end 
              
nvelo_si(npid,:)               = nmomentum_si(npid,:)/nmass_si(npid); 
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
%   % Update effective stress
for c=1:cellCount
    mpts = mspoints{c};
    dEFSP = zeros(2,2*length(mpts));
    dS = zeros(2*length(mpts),16);
    
    for sp = 1:length(mpts)
        spid = mpts(sp);
        L_sp = zeros(2,2);        
        
        for j=1:NODES(spid)
          if N{spid}(j)==0
         continue
          end
            npid = CONNECT{spid}(j);
            L_sp = L_sp + (nvelo_si(npid,:)'*dN{spid}(:,j)');
        end
              
        dESP = (L_sp + L_sp')/2*dt; 
        dEFSP(:,sp*2-1:sp*2) = dESP;  
        
        F_sp{spid} = (eye(2,2)+L_sp*dt)*F_sp{spid};
        
        r1_sp(spid,:) = (F_sp{spid} * r10_sp(spid,:)')';
        r2_sp(spid,:) = (F_sp{spid} * r20_sp(spid,:)')';        
                
        J = det(F_sp{spid});
        V_sp(spid)=V_spo(spid)*J;
         
        SSP = Lambda*log(J)/J*eye(2,2) + Mu/J*(F_sp{spid}*F_sp{spid}'-eye(2,2));
        s_sp(spid,1) = SSP(1,1);
        s_sp(spid,2) = SSP(2,2);
        s_sp(spid,3) = SSP(1,2);
        
%         de_sp(spid,1) = dEFSP(1,sp*2-1);
%         de_sp(spid,2) = dEFSP(2,sp*2);
%         de_sp(spid,3) = dEFSP(2,sp*2-1);
%         
%         de_sp(spid,1) = dESP(1,1);
%         de_sp(spid,2) = dESP(2,2);
%         de_sp(spid,3) = dESP(2,1);
%         
%         ds_sp(spid,:) = D * de_sp(spid,:)';
%         s_sp(spid,:)  = s_sp(spid,:) + ds_sp(spid,:);
    end
end

%% Update time and step
 t = t+dt;
 
 % Store old nodal acceleration
 nacc_old_si = nacc_si;
 nforce_old_si = nforce_si;
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