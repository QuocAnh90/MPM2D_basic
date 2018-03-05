% close all;
clear all;
tic;
addpath('F:/Git/MPM2D/SubFunctions/');
% Unit
% kN - seconds - m
%% Material porperties
E                       = 0.04;%10000          ;                  % Young modulus of solid
psp                     = 1.0;%4            ;                  % solid density
nu                      = 0;%.41           ;                  % Poison ratio
g                       = 0.0            ;                  % gravity acceleration
Lambda                  = E*nu/(1+nu)/(1-2*nu);
Mu                      = E/2/(1+nu);
Length                  = 4.0            ;   
%% Structured Grid input
NN(1)                 = 501;%33;                               % number of nodes in X direction
NN(2)                 = 4;%13;                               % number of nodes in Y direction
le(1)                 = 0.008;%0.5;                              % size of element in X direction
le(2)                 = le(1);                          % size of element in Y direction

%% Grid generation
[LOC,LOCC,cellCount,nodeCount] = Grid_Generation(NN,le);
% nodeCount: total number of nodes
% cellCount: total number of elements
% LOC(n,1:2): localtion coordinate of node "n" in x(1) and y(2) direction
% LOCC(e,1:2): localtion coordinate of element centroid "e" in x(1) and y(2) direction

%% Time
ftime                   = 5;%10;                               % Final time
dt                      = 0.01;%0.001;                            % Time step
ndt                     = round(ftime/dt) +1;
t                       = 0;

%% Boundary nodes
[nfbcx,nfbcy,fbcx,fbcy]=Compute_Boundary_Nodes(nodeCount,LOC);
% nfbcx: number of boundary nodes in X direction
% nfbcy: number of boundary nodes in Y direction
% fbcx: index of all boundary nodes in X direction
% fbcy: index of all boundary nodes in Y direction

%% Particle generation
% [spCount,x_sp,lp]=Particle_Ball(le);
% x_sp: Vector, position of MPs
% spCount: total number of MPs

% particle_per_cell       = 9;
% spCount                 = 16*particle_per_cell;                             % size of particle in Y direction
% lp(1)                 = le(1)/sqrt(particle_per_cell);                  % size of particle in X direction
% lp(2)                 = lp(1);  
% x_sp                    = zeros(spCount,2);
% sp=1;
% while sp<spCount+0.0001
%     for i=1:24
%         for j=1:6
%             x_sp(sp,1:2) = [2*le(1)+0.5*lp(1)+(i-1)*lp(1,1) 8*le(2)+0.5*lp(2)+(j-1)*lp(2)];
%             sp=sp+1;
%         end
%     end
% end

fid=fopen('vibrationx.txt');
particle_per_cell       = 4;
spCount                 = 498*particle_per_cell;
lp(1,1)                 = 0.004;                                 % size of particle in X direction
lp(1,2)                 = lp(1,1);                              % size of particle in Y direction
x_sp                    = zeros(spCount,1);
d_sp                    = zeros(spCount,2);

for sp=1:spCount
    x_sp(sp,1:2) = fscanf(fid,'%f',2);
end

%% Plot initial condition
initial_figure = Plot_Initial(x_sp,LOC,le);

%% Particle variables
dparticle               = lp(1)*lp(2);                          % area of particle domain
x_spo                   = x_sp;                                 % initial position
d_sp                    = zeros(spCount,2);                     % displacement
V_sp                    = dparticle * ones(spCount,1);          % volumn
V_spo                   = V_sp;                                 % initial volumn
m_sp                    = psp * V_sp;                           % mass
p_sp                    = psp * ones(spCount,1); 
b_sp                    = [0 0];                                % body force
s_sp                    = zeros(spCount,3);                     % Stress tensor
ds_sp                   = zeros(spCount,3);                     % Stress increment
v_ssp                   = zeros(spCount,2);                     % velocty
e_sp                    = zeros(spCount,3);                     % Strain tensor
de_sp                   = zeros(spCount,3);                     % Strain increment
ptraction_sp            = zeros(spCount,2);                     % traction
F_sp                    = cell(spCount,1);                      % Gradient deformation

%% Node variables
nmass_si                = zeros(nodeCount,1);                   % Nodal Mass
nmomentum_si            = zeros(nodeCount,2);                   % Nodal Momentum
niforce_si              = zeros(nodeCount,2);                   % Nodal Internal force
neforce_si              = zeros(nodeCount,2);                   % Nodal External force
nvelo_si                = zeros(nodeCount,2);                   % Nodal Velocity
traction_si             = zeros(nodeCount,2);                   % Nodal Traction

%% Initial condition
% Gradient deformation
for sp = 1:spCount
    F_sp{sp} = [1 0; 0 1];
end

% Velocity condition
% for sp=1:spCount
%     v_ssp(sp,:) = [5 0];
% end

 for sp=1: spCount
    e_sp(sp,1)               = 0.12 * (Length/2-x_sp(sp,1)) * exp(-60*(x_sp(sp,1)-Length/2).^2);
    s_sp(sp,1)               = e_sp(sp,1) * E;
    
    F_sp{sp}(1,1)               = F_sp{sp}(1,1) + e_sp(sp,1);
    end   

    ep_exact          = zeros(spCount,3);
    sp_exact        = zeros(spCount,3);
    
%     for sp=1:spCount
%     ep_exact(sp,1)         = 0.0005 * ((-120*x_sp(sp,1)+120*sqrt(E)*ftime+60*Length)*exp(-60*x_sp(sp,1).^2+sqrt(E)*ftime*(120*x_sp(sp,1)-60*Length)+60*x_sp(sp,1)*Length-60*E*ftime*ftime-15*Length*Length)+(-120*x_sp(sp,1)-120*sqrt(E)*ftime+60*Length)*exp(-60*x_sp(sp,1).^2-sqrt(E)*ftime*(120*x_sp(sp,1)-60*Length)+60*x_sp(sp,1)*Length-60*E*ftime*ftime-15*Length*Length));
%     sp_exact(sp,1)         = ep_exact(sp,1) * E;
%     end
    
%% start the algorithm
% video
timestep = 100;     % number of frame to save
r=timestep/20;      % number of frame per second video ~200s

writerObj2           = VideoWriter('MPM.avi');
writerObj2.FrameRate = r;    % number of frame per second
open(writerObj2);

    for tt = 1:timestep
    ft              = ftime/timestep*tt;
%     ft=ftime;
 while t<ft+0.0000000001      
     t
     
%% Reset value of node
nmass_si(:)                = 0;
nmomentum_si(:)            = 0;  
niforce_si(:)              = 0;
neforce_si(:)              = 0;
nvelo_si(:)                = 0;
traction_si(:)             = 0;

%% Store particles into cell
[N,dG,CONNECT,spElems,mspoints] = Compute_Interpolator(spCount,cellCount,x_sp,le,NN,LOC);
% N: Shape function
% dG: Gradient of shape function
% spElems(p): element index where "p" locates
% mspoints(e): all particles indexes where locate in the cell "e"

%% Mapping from particle to nodes using taylor least square (TLS)
% Compute Gauss point
GaussCount = 4;             % Number of Gauss Point per element
[x_gauss,y_gauss,omega_gauss,active_elements] = Locate_Gauss_Point(GaussCount,spElems,LOCC,le);

% Compute Conserve Integral
[Conserve_mass,Conserve_momentum] = Conserve_Integral(le,active_elements,mspoints,m_sp,v_ssp);

% Reconstruct function using TLS
[Reconstruction_density,Reconstruction_momentumx,Reconstruction_momentumy,Reconstruction_stressx,Reconstruction_stressy,Reconstruction_shear] = Reconstruction(active_elements,mspoints,p_sp,v_ssp,s_sp,GaussCount,x_gauss,y_gauss,x_sp,le,LOCC,Conserve_mass,Conserve_momentum);

% Gauss integration
[nmass_si,nmomentum_si,niforce_si] = Gauss_Integration(Reconstruction_density,Reconstruction_momentumx,Reconstruction_momentumy,Reconstruction_stressx,Reconstruction_stressy,Reconstruction_shear,CONNECT,LOC,le,GaussCount,x_gauss,y_gauss,omega_gauss,spElems,active_elements,nmass_si,nmomentum_si,niforce_si);

 %% Mapping from particle to nodes using MPM
% [nmass_si,nmomentum_si,niforce_si,neforce_si,traction_si]=Interpolate_Particle_To_Grid(CONNECT,le,N,dG,spCount,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,nmass_si,nmomentum_si,niforce_si,neforce_si,traction_si);
% nmass_si: nodal mass
% nmomentum_si: nodal momentum
% niforce_si: nodal internal force
% neforce_si: nodal external force
% traction_si: nodal traction

%% Update momentum
% Update force and momentum
 nforce_si      = niforce_si + neforce_si + traction_si;
 nmomentum_si   = nmomentum_si + nforce_si*dt;
 
% Boundary condition
[nforce_si]     = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nforce_si); % Boundary condition for nodal force
[nmomentum_si]  = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nmomentum_si); % Boundary condition for nodal force

%% Update solid particle velocity and position
[v_ssp,x_sp,d_sp] = Update_Particle_Position(dt,CONNECT,N,spCount,nmass_si,nforce_si,nmomentum_si,x_spo,v_ssp,x_sp);
% velocity particle: v_ssp
% position particle: x_sp
% displacement particle: d_sp
 
%% Mapping nodal velocity back to node
[nvelo_si] = Interpolate_velocity_back(spCount,CONNECT,nmomentum_si,m_sp,v_ssp,N,nmass_si);
% Boundary condition
[nvelo_si] = Boundary_Dirichlet(nfbcx,nfbcy,fbcx,fbcy,nvelo_si); % Boundary condition for nodal force

%% Update effective stress
[F_sp,V_sp,s_sp,p_sp] = Update_Stress(dt,cellCount,mspoints,CONNECT,nvelo_si,dG,Lambda,Mu,F_sp,V_spo,m_sp,s_sp,p_sp,V_sp);
% F_sp: gradient deformation
% V_sp: volumn
% s_sp: Cauchy stress
% p_sp: density

 % Update time and step 
 t = t+dt;
 end

 %% Plot the result
%     StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,e_sp,v_ssp,spCount);

    for sp=1:spCount
    ep_exact(sp,1)         = 0.0005 * ((-120*x_sp(sp,1)+120*sqrt(E)*t+60*Length)*exp(-60*x_sp(sp,1).^2+sqrt(E)*t*(120*x_sp(sp,1)-60*Length)+60*x_sp(sp,1)*Length-60*E*t*t-15*Length*Length)+(-120*x_sp(sp,1)-120*sqrt(E)*t+60*Length)*exp(-60*x_sp(sp,1).^2-sqrt(E)*t*(120*x_sp(sp,1)-60*Length)+60*x_sp(sp,1)*Length-60*E*t*t-15*Length*Length));
    sp_exact(sp,1)         = ep_exact(sp,1) * E;
    end

    StressProfile1=figure;
    set(StressProfile1, 'visible','off');
    % plot(x_sp(:,1),e_sp(:,1),'.',x_sp(:,1),ep_exact(:,1));
    plot(x_sp(:,1),s_sp(:,1),'.',x_sp(:,1),sp_exact(:,1));
    axis([2,max(LOC(:,1)),-0.0002,0.0002]);
    legend('2D-MPM-TLS','analytical solution')
    title('Strain plot of TLS-MPM solution in 2D'); % title;
    ylabel('Strain'); % label for y axis
    xlabel('Length (m)'); % label for x axis
    
    frame2 = getframe(StressProfile1);
    writeVideo(writerObj2,frame2);
    end
    
    close(writerObj2);