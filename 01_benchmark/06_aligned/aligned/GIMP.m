close all;
clear;
tic;

% Unit
% kN - seconds - m
%% Material porperties
E                       = 10000000     ;                  % Young modulus of solid
psp                     = 1000           ;                  % solid density
nu                      = 0.3           ;                  % Poison ratio
g                       = 0.0            ;                  % gravity acceleration

C                       = sqrt(E/psp);
D = [E/(1-nu^2) E*nu/(1-nu^2) 0; E*nu/(1-nu^2) E/(1-nu^2) 0; 0 0 E/(1-nu)];

Lambda                  = E*nu/(1+nu)/(1-2*nu);
Mu                      = E/2/(1+nu);

%% Grid generation
resolution              = 64;

NN(1,1)                 = 1.5/(1/resolution)+1;               % number of nodes in X direction
NN(1,2)                 = NN(1,1);                               % number of nodes in Y direction
cellCount               = (NN(1)-1)*(NN(2)-1);              % number of elements
nodeCount               = NN(1)*NN(2);                      % number of nodes
LOC                     = zeros(NN(1)*NN(2),2);             % zero matrix of all nodal coordinate
le(1,1)                 = 1/resolution;                              % size of element in X direction
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

%% Time
ftime                   = 0.02;
dt                      = 0.4*le(1,1)/C;
ndt                     = round(ftime/dt) +1;
t                       = 0;

% time integration
p_b                     = 1.0;
a_m                     = (2*p_b-1)/(1+p_b);
beta                    = (5-3*p_b)/(1+p_b).^2/(2-p_b);
gamma                   = 3/2-a_m;

%% Boundary condition
nfbc                    = 0            ;                  % number of fixed nodes
nfbn                    = 0            ;                  % number of traction nodes
fbc =[];

nfbcx = 0;
nfbcy = 0;
fbcx =[];
fbcy =[];
for n=1:nodeCount
    if LOC(n,1)<=0.25 || LOC(n,1)>=1.25
        nfbcx = nfbcx + 1;
        fbcx = [fbcx n];
    end
end

for n=1:nodeCount
    if LOC(n,2)<=0.25 || LOC(n,2)>=1.25
        nfbcy = nfbcy + 1;
        fbcy = [fbcy n];
    end
end

%% Node Variables
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

%% Particle generation
particle_per_cell       = 4;
spCount                 = resolution*resolution*particle_per_cell;
lp(1,1)                 = le(1,1)/2;                                 % size of particle in X direction
lp(1,2)                 = lp(1,1);                              % size of particle in Y direction
x_sp                    = zeros(spCount,2);
d_sp                    = zeros(spCount,2);

sp=1;
while sp<spCount+0.0001
    for i=1:resolution*2
        for j=1:resolution*2
            x_sp(sp,1:2)= [0.25+0.5*lp(1,1)+(j-1)*lp(1,1) 0.25+0.5*lp(1,2)+(i-1)*lp(1,2)];
            sp=sp+1;
        end
    end  
end

dparticle               = lp(1)*lp(2);                          % area of particle domain
x_spo                   = x_sp;                                 % initial position
b_sp                    = zeros(spCount,2);                                % body force
s_sp                    = zeros(spCount,3);                     % Stress tensor
ds_sp                   = zeros(spCount,3);                     % Stress increment
v_ssp                   = zeros(spCount,2);                     % velocty
e_sp                    = zeros(spCount,3);                     % Strain tensor
de_sp                   = zeros(spCount,3);                     % Strain increment
% F_sp                    = zeros(spCount,3);
ptraction_sp            = zeros(spCount,2);
alpha_sp                = zeros(spCount,1);
r1_sp                   = zeros(spCount,2);
r2_sp                   = zeros(spCount,2);
 x_corner1              = zeros(spCount,2);
 x_corner2              = zeros(spCount,2);
 x_corner3              = zeros(spCount,2);
 x_corner4              = zeros(spCount,2);
 
 CONNECT                = zeros(spCount,16);
 N                      = zeros(spCount,16);
 dN                     = zeros(spCount,32);
 dG                     = zeros(16*spCount,2);
 spElems                = zeros(spCount,1); 
 mspoints               = cell(cellCount,1);
 F_sp                   = cell(spCount,1);
  
deviation_X               = zeros(spCount,2);
deviation_V               = zeros(spCount,2);
deviation_S               = zeros(spCount,3);
dis_ana                 = zeros(spCount,2);
vel_ana                 = zeros(spCount,2);
stress_ana                 = zeros(spCount,3);

ux                      = zeros(spCount,1);
uy                      = zeros(spCount,1);
vx                      = zeros(spCount,1);
vy                      = zeros(spCount,1);
sx                      = zeros(spCount,1);
sy                      = zeros(spCount,1);
Fxx                     = zeros(spCount,1);
Fyy                     = zeros(spCount,1);

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

%% Initial condition
A = 0.1;
% A = 0.001;
sx1=[];
sy1=[];
sx2=[];
sy2=[];
ux1=[];
uy1=[];
time=[];
error_X=[];
error_V=[];
error_S=[];
perror_X=0;
perror_V=0;
perror_S=0;
dx1=[];
dy1=[];
Fxx1=[];
Fxx2=[];
Fyy1=[];
Fyy2=[];
bxx=[];
byy=[];

for sp=1:spCount
    v_ssp(sp,1) = A*C*pi*sin(2*pi*(x_spo(sp,1)-0.25))*cos(C*pi*t);
    v_ssp(sp,2) = A*C*pi*sin(2*pi*(x_spo(sp,2)-0.25))*cos(C*pi*t+pi); 
end

% %% Visualization
%  for sp=1:spCount
%  x_corner1(sp,:) = x_sp(sp,:) - r1_sp(sp,:) - r2_sp(sp,:);      % Position of corner 1
%  x_corner2(sp,:) = x_sp(sp,:) + r1_sp(sp,:) - r2_sp(sp,:);
%  x_corner3(sp,:) = x_sp(sp,:) + r1_sp(sp,:) + r2_sp(sp,:);
%  x_corner4(sp,:) = x_sp(sp,:) - r1_sp(sp,:) + r2_sp(sp,:);
%  end
%  
%  figure
% plot(x_sp(:,1),x_sp(:,2),'.');
% title('initial condition');
% grid on
% hold on
% for sp=1:spCount
%     x_cor = [x_corner1(sp,1) x_corner2(sp,1) x_corner3(sp,1) x_corner4(sp,1) x_corner1(sp,1)];
%     y_cor = [x_corner1(sp,2) x_corner2(sp,2) x_corner3(sp,2) x_corner4(sp,2) x_corner1(sp,2)];
%     plot(x_cor,y_cor,'r')
% end
% axis([0,1.5,0,1.5]);
% set(gca,'xtick',(0:le(1,1):LOCX(NN(1,1))));
% set(gca,'ytick',(0:le(1,2):LOCY(NN(1,2))));
% % set(gca,'XTickLabel',[]);
% % set(gca,'YTickLabel',[]);

% %% start the algorithm
% % video
% timestep = 100;    % number of frame to save
% r=timestep/20;      % number of frame per second video ~200s
% 
% writerObj2           = VideoWriter('gimp.avi');
% writerObj2.FrameRate = r;    % number of frame per second
% open(writerObj2);
%     for tt = 1:timestep
%          ft              = ftime/timestep*tt;
         
  ft=ftime;

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

for sp=1:spCount
    ux(sp) = A*sin(2*pi*(x_spo(sp,1)-0.25))*sin(C*pi*t);
    uy(sp) = A*sin(2*pi*(x_spo(sp,2)-0.25))*sin(C*pi*t+pi);
    
    vx(sp) = A*C*pi*sin(2*pi*(x_spo(sp,1)-0.25))*cos(C*pi*t);
    vy(sp) = A*C*pi*sin(2*pi*(x_spo(sp,2)-0.25))*cos(C*pi*t+pi);
    
    Fxx(sp) = 1 + 2*A*pi*cos(2*pi*(x_spo(sp,1)-0.25))*sin(C*pi*t);
    Fyy(sp) = 1 + 2*A*pi*cos(2*pi*(x_spo(sp,2)-0.25))*sin(C*pi*t+pi);
    
    sx(sp) = (Lambda*log(Fxx(sp)*Fyy(sp))+Mu*(Fxx(sp)^2-1))/Fxx(sp)/Fyy(sp);
    sy(sp) = (Lambda*log(Fxx(sp)*Fyy(sp))+Mu*(Fyy(sp)^2-1))/Fxx(sp)/Fyy(sp);
    
%     SSP = Lambda*log(J)/J*eye(2,2) + Mu/J*(F_sp{spid}*F_sp{spid}'-eye(2,2));
    
    dis_ana(sp,:) =  [ux(sp) uy(sp)];
    vel_ana(sp,:) =  [vx(sp) vy(sp)];
    stress_ana(sp,:) =  [sx(sp) sy(sp) 0];
    
    b_sp(sp,1) = pi^2*ux(sp)*(4*Mu/psp-C^2-4*(Lambda*(log(Fxx(sp)*Fyy(sp))-1)-Mu)/psp/Fxx(sp)^2);
    b_sp(sp,2) = pi^2*uy(sp)*(4*Mu/psp-C^2-4*(Lambda*(log(Fxx(sp)*Fyy(sp))-1)-Mu)/psp/Fyy(sp)^2);   
end
Fxx1 = [Fxx1 Fxx(1)];
Fyy1 = [Fyy1 Fyy(2)];
bxx = [bxx b_sp(1,1)];
byy = [byy b_sp(1,2)];

%% Store particles into cell
 for sp=1:spCount
 x_corner1(sp,:) = x_sp(sp,:) - r1_sp(sp,:) - r2_sp(sp,:);      % Position of corner 1
 x_corner2(sp,:) = x_sp(sp,:) + r1_sp(sp,:) - r2_sp(sp,:);
 x_corner3(sp,:) = x_sp(sp,:) + r1_sp(sp,:) + r2_sp(sp,:);
 x_corner4(sp,:) = x_sp(sp,:) - r1_sp(sp,:) + r2_sp(sp,:);
 end
 
 for sp = 1:spCount
 spElems(sp) = ceil(x_sp(sp,1)/le(1))+(NN(1)-1)*(fix(x_sp(sp,2)/le(2)));   % compute vector store index elements                           
 
 CONNECT(sp,6)  = spElems(sp) + floor(spElems(sp)/(NN(1)-1));
 CONNECT(sp,2)  = CONNECT(sp,6)-NN(1);
 CONNECT(sp,1)  = CONNECT(sp,2)-1;
 CONNECT(sp,3)  = CONNECT(sp,2)+1;
 CONNECT(sp,4)  = CONNECT(sp,2)+2;
 CONNECT(sp,5)  = CONNECT(sp,6)-1;
 CONNECT(sp,7)  = CONNECT(sp,6)+1;
 CONNECT(sp,8)  = CONNECT(sp,6)+2;
 CONNECT(sp,10) = CONNECT(sp,6)+NN(1);
 CONNECT(sp,11) = CONNECT(sp,7)+NN(1);
 CONNECT(sp,9)  = CONNECT(sp,10)-1;
 CONNECT(sp,12) = CONNECT(sp,11)+1;
 CONNECT(sp,13) = CONNECT(sp,9)+NN(1);
 CONNECT(sp,14) = CONNECT(sp,10)+NN(1);
 CONNECT(sp,15) = CONNECT(sp,11)+NN(1);
 CONNECT(sp,16) = CONNECT(sp,12)+NN(1);
 
     for i=1:16
[N(sp,i),dN(sp,i),dN(sp,i+16)]=GIMPshape(x_sp(sp,1:2),LOC(CONNECT(sp,i),:),le(1,1),le(1,2),lp(1,1));
     end
 
 % Build matrix of gradient of shape function of 4 nodes of the cell which
 % contain the particles
 for i=1:16
dG((sp-1)*16+i,:) = [dN(sp,i) dN(sp,i+16)];
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
 
 for j=1:16
     npid                           = CONNECT(sp,j);
     
 % Mass
 nmass_si(npid)            = nmass_si(npid) + m_sp(sp)*N(sp,j);
 
 % Momentumlp(1,1)lp(1,1)lp(1,1)
 nmomentum_si(npid,:)      = nmomentum_si(npid,:) + m_sp(sp)*v_ssp(sp,:)*N(sp,j);
 
 % Internal forces
niforce_si(npid,:)         = niforce_si(npid,:) - (V_sp(sp)*SSP*dG((sp-1)*16+j,:)')';

 % External forces
neforce_si(npid,:)         = neforce_si(npid,:) +  b_sp(sp,:)*m_sp(sp)*N(sp,j);

% Traction
 traction_si(npid,:)       = traction_si(npid,:) + V_sp(sp)*ptraction_sp(sp,:)*N(sp,j)/le(1,1)/le(1,2);
 
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
end

for i=1:nfbcy
nforce_si(fbcy(i),2)     = 0;
nmomentum_si(fbcy(i),2)  = 0;
end

 %% Update solid particle velocity and position
 % Update position and velocity
 for sp = 1:spCount

     for j = 1:16
              if nmass_si(CONNECT(sp,j))==0
                continue
              end         
         v_ssp(sp,:)                      = v_ssp(sp,:) + dt * N(sp,j) * (nforce_si(CONNECT(sp,j),:)/nmass_si(CONNECT(sp,j)));
         x_sp(sp,:)                       = x_sp(sp,:) + nmomentum_si(CONNECT(sp,j),:)*N(sp,j)*dt/ nmass_si(CONNECT(sp,j));
         d_sp(sp,:)                       = d_sp(sp,:) + nmomentum_si(CONNECT(sp,j),:)*N(sp,j)*dt/ nmass_si(CONNECT(sp,j));
     end   
 end
 
 %% Mapping velocity back to node
 nmomentum_si(:) = 0;
 % Momentum
 for sp = 1:spCount
     for j = 1:16
         nmomentum_si(CONNECT(sp,j),:)  = nmomentum_si(CONNECT(sp,j),:) + m_sp(sp)*v_ssp(sp,:)*N(sp,j);
     end
 end
 
 % Velocity
  for sp = 1:spCount
 for j = 1:16
              if nmass_si(CONNECT(sp,j))==0
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

%   % Update effective stress
for c=1:cellCount
    mpts = mspoints{c};
    dEFSP = zeros(2,2*length(mpts));
    dS = zeros(2*length(mpts),16);
    
    for sp = 1:length(mpts)
        spid = mpts(sp);
        L_sp = zeros(2,2);
    
        for j=1:16
            L_sp = L_sp + (nvelo_si(CONNECT(spid,j),:)'*dG((spid-1)*16+j,:));
        end
        
%         dESP = (L_sp + L_sp')/2*dt; 
%         dEFSP(:,sp*2-1:sp*2) = dESP;

        F_sp{spid} = (eye(2,2)+L_sp*dt)*F_sp{spid};
        
        r1_sp(spid,:) = (F_sp{spid} * r10_sp(spid,:)')';
        r2_sp(spid,:) = (F_sp{spid} * r20_sp(spid,:)')';   
        
        J = det(F_sp{spid});
        V_sp(spid)=V_spo(spid)*J;
        
        SSP = Lambda*log(J)/J*eye(2,2) + Mu/J*(F_sp{spid}*F_sp{spid}'-eye(2,2));
        s_sp(spid,1) = SSP(1,1);
        s_sp(spid,2) = SSP(2,2);
        s_sp(spid,3) = SSP(1,2);
        
    end
       
%     for sp = 1:length(mpts)
%         spid = mpts(sp);
%         de_sp(spid,1) = dEFSP(1,sp*2-1);
%         de_sp(spid,2) = dEFSP(2,sp*2);
%         de_sp(spid,3) = dEFSP(2,sp*2-1);
%         
%         ds_sp(spid,:) = D * de_sp(spid,:)';
%         s_sp(spid,:)  = s_sp(spid,:) + ds_sp(spid,:);
%     end
end

%% Compute error
x=0;
v=0;
s=0;

for sp=1:spCount
    deviation_X(sp,:) = dis_ana(sp,:)-d_sp(sp,:);
    x=x+(norm(deviation_X(sp,:)))^2;
    
    deviation_V(sp,:) = vel_ana(sp,:)-v_ssp(sp,:);
    v=v+(norm(deviation_V(sp,:)))^2;
    
    deviation_S(sp,:) = stress_ana(sp,:)-s_sp(sp,:);
    s=s+(norm(deviation_S(sp,:)))^2;
end
    perror_X = perror_X + x;
    error_X = [error_X sqrt(x)/spCount];

    perror_V = perror_V + v;
    error_V = [error_V sqrt(v)/spCount];
    
    perror_S = perror_S + s;
    error_S = [error_S sqrt(s)/spCount];
    
    ux1 = [ux1 ux(1)];
    uy1 = [uy1 uy(1)];
    time = [time t];
    sx1 = [sx1 sx(1)];
    sy1 = [sy1 sy(1)];
    sx2 = [sx2 s_sp(1,1)];
    sy2 = [sy2 s_sp(1,2)];    
    
    dx1 = [dx1 d_sp(1,1)];
    dy1 = [dy1 d_sp(1,2)];
    
    Fxx2 = [Fxx2 F_sp{1}(1,1)];
    Fyy2 = [Fyy2 F_sp{1}(2,2)];
    
 % Update time and step 
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

 for sp=1:spCount
 x_corner1(sp,:) = x_sp(sp,:) - r1_sp(sp,:) - r2_sp(sp,:);      % Position of corner 1
 x_corner2(sp,:) = x_sp(sp,:) + r1_sp(sp,:) - r2_sp(sp,:);
 x_corner3(sp,:) = x_sp(sp,:) + r1_sp(sp,:) + r2_sp(sp,:);
 x_corner4(sp,:) = x_sp(sp,:) - r1_sp(sp,:) + r2_sp(sp,:);
 end
 
StressProfile1=figure;
set(StressProfile1, 'visible','off');
sz = 2;
color = strain;
% color = s_sp(:,3);
scatter(x_sp(:,1),x_sp(:,2),sz,color,'filled');
hold on
for sp=1:spCount
    x_cor = [x_corner1(sp,1) x_corner2(sp,1) x_corner3(sp,1) x_corner4(sp,1) x_corner1(sp,1)];
    y_cor = [x_corner1(sp,2) x_corner2(sp,2) x_corner3(sp,2) x_corner4(sp,2) x_corner1(sp,2)];
    plot(x_cor,y_cor,'r')
end
grid on
axis([0,LOCX(NN(1,1)),0,LOCY(NN(1,2))]);
set(gca,'xtick',(0:0.5:LOCX(NN(1,1))));
set(gca,'ytick',(0:0.5:LOCY(NN(1,2))));
colormap(jet(256))
% set(h, 'ytick', [0:0.2:1.2]);
caxis([0 0.25]);

figure
plot(time, ux1, time, dx1);
title('ux');

figure
plot(time, uy1, time, dy1);
title('uy');

figure 
plot(time,Fxx1,time,Fxx2);
title('Fx');

figure
plot(time,Fyy1,time,Fyy2);
title('Fy');

% figure
% plot(time,bxx);
% title('bx');
% 
% figure
% plot(time,byy);
% title('by');

%     frame2 = getframe(StressProfile1);
%     writeVideo(writerObj2,frame2);
%      end
%     close(writerObj2);
%     
    perror_X = sqrt(perror_X)/ndt/spCount;
    error_X=error_X';
    perror_V = sqrt(perror_V)/ndt/spCount;
    error_V=error_V';    
    perror_S = sqrt(perror_S)/ndt/spCount;
    error_S=error_S';        
    
    time=time';