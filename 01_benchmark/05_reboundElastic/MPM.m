% close all;
clear all;
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

%% Boundary condition
nfbc                    = 82             ;                  % number of fixed nodes
nfbn                    = 0              ;                  % number of traction nodes
fbc = [];

for i=1:31
fbc                     = [fbc 35+i-1]         ;                  % index of fixed nodes
end

for i=1:31
fbc                     = [fbc 365+i-1]         ;                  % index of fixed nodes
end

for i=1:10
fbc                     = [fbc 35+33*i]         ;                  % index of fixed nodes
end

for i=1:10
fbc                     = [fbc 65+33*i]         ;                  % index of fixed nodes
end

%% Time
ftime                   = 10;
dt                      = 0.001;
ndt                     = round(ftime/dt) +1;
t                       = 0;
n                       = 1;
damps                   = 0.0;

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

LOCX                    = [0:NN(1)-1]'*le(1);               % Location of all nodes in X direction
LOCY                    = [0:NN(2)-1]'*le(2);               % Location of all nodes in Y direction

LOCCX                   = [0:(NN(1)-1)-1]'*le(1)+le(1)/2;   % Location of cells in X direction
LOCCY                   = [0:(NN(2)-1)-1]'*le(2)+le(2)/2;   % Location of cells in X direction

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
set(gca,'xtick',[0:le(1,1):LOCX(NN(1,1))]);
set(gca,'ytick',[0:le(1,2):LOCY(NN(1,2))]);
% set(gca,'XTickLabel',[]);
% set(gca,'YTickLabel',[]);

%% Initial condition
% Particle
dparticle               = lp(1)*lp(2);                          % area of particle domain
x_spo                   = x_sp;                                 % initial position
V_sp                    = dparticle * ones(spCount,1);          % volumn
m_sp                    = psp * V_sp;                           % mass
b_sp                    = [0 0];                                % body force
s_sp                    = zeros(spCount,3);                     % Stress tensor
ds_sp                   = zeros(spCount,3);                     % Stress increment
v_ssp                   = zeros(spCount,2);                     % velocty
e_sp                    = zeros(spCount,3);                     % Strain tensor
de_sp                   = zeros(spCount,3);                     % Strain increment
F_sp                    = zeros(spCount,3);
ptraction_sp            = zeros(spCount,2);
alpha_sp                = zeros(spCount,1);

for sp=1:spCount
    v_ssp(sp,:) = [5 0];
end

% Node
nmass_si                = zeros(nodeCount,1);
nmomentum_si            = zeros(nodeCount,2);
nvolume_si              = zeros(nodeCount,1);
niforce_si              = zeros(nodeCount,2);
neforce_si              = zeros(nodeCount,2);
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

writerObj2           = VideoWriter('MPM.avi');
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
nacc_si(:)                 = 0;
nacc_middle_si(:)          = 0;
traction_si(:)             = 0;

%% Store particles into cell
 spElems = zeros(spCount,1);                            % vectors store no of elements for each particle
 xc = zeros(spCount,2);                                 % particle position - centroid of element
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
dG((sp-1)*4+[1:4],:) = [dN(sp,1) dN(sp,5);dN(sp,2) dN(sp,6);dN(sp,3) dN(sp,7);dN(sp,4) dN(sp,8)];
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
     spid                           = CONNECT(sp,j);
     
 % Mass
 nmass_si(CONNECT(sp,j))            = nmass_si(CONNECT(sp,j)) + m_sp(sp)*N(sp,j);
 
 % Momentum
 nmomentum_si(CONNECT(sp,j),:)      = nmomentum_si(CONNECT(sp,j),:) + m_sp(sp)*v_ssp(sp,:)*N(sp,j);
 
 % Internal forces
niforce_si(CONNECT(sp,j),:)         = niforce_si(CONNECT(sp,j),:) - (V_sp(sp)*SSP*dG((sp-1)*4+j,:)')';

 % External forces
neforce_si(CONNECT(sp,j),:)         = neforce_si(CONNECT(sp,j),:) + b_sp*m_sp(sp)*N(sp,j);

% Traction
 traction_si(CONNECT(sp,j),:)       = traction_si(CONNECT(sp,j),:) + V_sp(sp)*ptraction_sp(sp,:)*N(sp,j)/le(1,1)/le(1,2);
 
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
 nacc_si(n,:)        = (nacc_middle_si(n,:) - nacc_old_si(n,:) * a_m)/(1 - a_m);
 end
 
 % Boundary condition
for i=1:nfbc
nforce_si(fbc(i),:)     = 0;
nmomentum_si(fbc(i),:)  = 0;
end
 
 % Update solid particle velocity and position
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
for i=1:nfbc
nvelo_si(fbc(i),:)     = 0;
end
  
%   % Update effective stress
for c=1:cellCount
    mpts = mspoints{c};
    dEFSP = zeros(2,2*length(mpts));
    dS = zeros(2*length(mpts),4);
    
    for sp = 1:length(mpts)
        spid = mpts(sp);
        L_sp = zeros(2,2);
    
        for j=1:4
            L_sp = L_sp + [nvelo_si(CONNECT(spid,j),:)'*dG((spid-1)*4+j,:)];
        end
        
        dESP = (L_sp + L_sp')/2*dt; 
        dEFSP(:,sp*2-1:sp*2) = dESP;
        dS(2*sp-1:2*sp,:) = dG((spid-1)*4+1:(spid-1)*4+4,:)';
    end
    
    if cfilter==1
        
    % QR method
    if isempty(mpts) ==0
        [Q1,R]=qr(dS,0);
        
        r = rank(R);
%         if r~=3
%             break
%         end
        dEFSP = (dEFSP*Q1(:,1:r))*Q1(:,1:r)';
    end    
    end
    
    for sp = 1:length(mpts)
        spid = mpts(sp);
        de_sp(spid,1) = dEFSP(1,sp*2-1);
        de_sp(spid,2) = dEFSP(2,sp*2);
        de_sp(spid,3) = dEFSP(2,sp*2-1);
        
        ds_sp(spid,:) = D * de_sp(spid,:)';
        e_sp(spid,:)  =  e_sp(spid,:)+de_sp(spid,:);
        s_sp(spid,:)  = s_sp(spid,:) + ds_sp(spid,:);
    end
end

 % Update time and step 
 n = n+1;
 t = t+dt;
 
 % Store old nodal acceleration
 for n=1:nodeCount
 nacc_old_si(n,:) = nacc_si(n,:);
 end
 
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
set(gca,'xtick',[0:0.5:LOCX(NN(1,1))]);
set(gca,'ytick',[0:0.5:LOCY(NN(1,2))]);
h=colorbar;
colormap(jet(256))

% set(h, 'ytick', [0:0.2:1.2]);
caxis([0 0.25]);

    frame2 = getframe(StressProfile1);
    writeVideo(writerObj2,frame2);
 
    end

    close(writerObj2);