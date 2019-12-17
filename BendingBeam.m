close all;
clear all;
tic;
addpath('SubFunctions/Cores');
addpath('SubFunctions/Constitutive_Models');
addpath('SubFunctions/CPDI_solvers');
addpath('SubFunctions/MPM_solvers');
addpath('SubFunctions/TLS_MPM_solvers');
addpath('SubFunctions/MLS_MPM_solvers');
addpath('SubFunctions/TLS_CPDI_solvers');
addpath('SubFunctions/TLS');
addpath('Particle_generation_input');

%  test
% test again

% Unit
% Newton - seconds - metre

%% Please select the versions of MPM!!!!!!!!!!!!!!!!!
% Original MPM: 'MPM'
% CPDI: 'CPDI'

% Please remember to change the output video file!!
Version = 'MPM';

%% Constutitive model
% Linear_Elastic
% Neo_Hookean_Elastic
CModel = 'Neo_Hookean_Elastic';

%% Material porperties
E                       = 1000000          ;                  % Young modulus of solid
psp                     = 1050.0            ;                  % solid density
nu                      = 0.3           ;                  % Poison ratio
g                       = 10.0            ;                  % gravity acceleration

CModel_parameter = [E,nu];

%% Structured Grid input
NN(1)                 = 13;%33;                               % number of nodes in X direction
NN(2)                 = 13;%13;                               % number of nodes in Y direction
le(1)                 = 0.5;%0.5;                              % size of element in X direction
le(2)                 = le(1);                          % size of element in Y direction

%% Grid generation
[LOC,LOCC,cellCount,nodeCount] = Grid_Generation(NN,le);
% nodeCount: total number of nodes
% cellCount: total number of elements
% LOC(n,1:2): localtion coordinate of node "n" in x(1) and y(2) direction
% LOCC(e,1:2): localtion coordinate of element centroid "e" in x(1) and y(2) direction

%% Time
c                       = sqrt(E/psp);
ftime                   = 5;
dt                      = 0.001;
ndt                     = round(ftime/dt) +1;
t                       = 0;

%% Boundary nodes
% Boundary coordination
x_min = 1;
x_max = 6;
y_min = 0;
y_max = 6;
[nfbcx,nfbcy,fbcx,fbcy]=Compute_Boundary_Nodes(nodeCount,LOC,x_max,x_min,y_max,y_min);
% nfbcx: number of boundary nodes in X direction
% nfbcy: number of boundary nodes in Y direction
% fbcx: index of all boundary nodes in X direction
% fbcy: index of all boundary nodes in Y direction

for n=1:nodeCount
    if LOC(n,1)<=x_min 
        nfbcy = nfbcy+1;
        fbcy = [fbcy n];
    end
end

%% Particle generation
particle_per_cell       = 9;
spCount                 = 16*particle_per_cell;
lp(1)                   = le(1)/sqrt(particle_per_cell);                                 % size of particle in X direction
lp(2)                   = lp(1);                              % size of particle in Y direction
x_sp                    = zeros(spCount,1);
d_sp                    = zeros(spCount,2);

sp=1;
while sp<spCount+0.0001
    for i=1:6
        for j=1:24
            x_sp(sp,1:2)= [2*le(1,1)+0.5*lp(1,1)+(j-1)*lp(1,1) 8*le(1,2)+0.5*lp(1,2)+(i-1)*lp(1,2)];
            sp=sp+1;
        end
    end
end
% x_sp: Vector, position of MPs
% spCount: total number of MPs

%% Plot initial condition
initial_figure = Plot_Initial(x_sp,LOC,le);

%% Particle variables
dparticle               = lp(1)*lp(2);                          % area of particle domain
x_spo                   = x_sp;                                 % initial position
p_sp                    = psp * ones(spCount,1);                % Density
b_sp                    = [zeros(spCount,1) -g*ones(spCount,1)];% body force

d_sp                    = zeros(spCount,2);                     % displacement
s_sp                    = zeros(spCount,3);                     % Stress tensor
ds_sp                   = zeros(spCount,3);                     % Stress increment
v_ssp                   = zeros(spCount,2);                     % velocty
e_sp                    = zeros(spCount,3);                     % Strain tensor
de_sp                   = zeros(spCount,3);                     % Strain increment
ptraction_sp            = zeros(spCount,2);                     % traction
F_sp                    = cell(spCount,1);                      % Gradient deformation
r1_sp                   = zeros(spCount,2);
r2_sp                   = zeros(spCount,2);

%% Initial condition
% Gradient deformation
for sp = 1:spCount
    r1_sp(sp,:) = [lp(1,1)/2 0];
    r2_sp(sp,:) = [0 lp(1,2)/2];
    F_sp{sp} = [1 0; 0 1];
end
r10_sp = r1_sp;
r20_sp = r2_sp;
V_sp                    = zeros(spCount,1);
for sp=1:spCount
V_sp(sp)                = 4*abs(r1_sp(sp,1)*r2_sp(sp,2)-r1_sp(sp,2)*r2_sp(sp,1)); 
end
V_spo                   = V_sp;
m_sp                    = psp * V_sp;                           % mass

% Traction
for sp=1:spCount
    if x_sp(sp,1)>1.01
        ptraction_sp(sp,1) = 1;
    end
end

%% start the algorithm
% video
timestep = 100;     % number of frame to save
r=timestep/20;      % number of frame per second video ~200s

writerObj2           = VideoWriter('BendingBeam.avi');
writerObj2.FrameRate = r;    % number of frame per second
open(writerObj2);

    for tt = 1:timestep
    ft              = ftime/timestep*tt;
%     ft=ftime;
 while t<ft+0.0000000001      
     t
     
     switch Version
         case 'MPM'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = MPM_solver(CModel,CModel_parameter,...
    nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
    nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt);

         case 'CPDI'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,r1_sp,r2_sp] = CPDI_solver(CModel,CModel_parameter,...
            nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
            nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp);
      end
     
     
 % Update time and step 
 t = t+dt;
 end

 %% Plot the result
StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,e_sp,v_ssp,spCount,r1_sp,r2_sp);

    frame2 = getframe(StressProfile1);
    writeVideo(writerObj2,frame2);
    end
    
    close(writerObj2);
  