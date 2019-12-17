% close all;
clear all;
tic;
addpath('SubFunctions/Cores');
addpath('SubFunctions/Constitutive_Models');
addpath('SubFunctions/CPDI_solvers');
addpath('SubFunctions/MPM_solvers');
addpath('SubFunctions/MLS_MPM_solvers');
addpath('SubFunctions/TLS_MPM_solvers');
addpath('SubFunctions/TLS_CPDI_solvers');
addpath('SubFunctions/KIMP_solvers');
addpath('SubFunctions/TLS');
addpath('Particle_generation_input');

% Unit
% Newton - seconds - metre

%% Please select the versions of MPM!!!!!!!!!!!!!!!!!
% Original MPM: 'MPM'
% B spline MPM: 'BSMPM'
% CPDI: 'CPDI'
% Moving least square MPM: 'MLS_MPM'

% Please remember to change the output video file!!
Version = 'MLS_MPM';

%% Constutitive model
% Linear_Elastic
% Neo_Hookean_Elastic
CModel = 'Neo_Hookean_Elastic';

% Amplitude of displacement
A = 0.0001;
particle_index = 2;

%% Material porperties
E                       = 10000000          ;                % Young modulus of solid
psp                     = 1000.0            ;               % solid density
nu                      = 0.3           ;                   % Poison ratio
g                       = 0.0            ;                  % gravity acceleration
Lambda                  = E*nu/(1+nu)/(1-2*nu);
Mu                      = E/2/(1+nu);
CModel_parameter = [E,nu];

%% Structured Grid input
resolution              = 5;
NN(1)                   = resolution + 5;                               % number of nodes in X direction
NN(2)                   = 4;                               % number of nodes in Y direction
le(1)                   = 1/resolution;                              % size of element in X direction
le(2)                   = 1/resolution;                            % size of element in Y direction

%% Grid generation
[LOC,LOCC,cellCount,nodeCount] = Grid_Generation(NN,le);
% nodeCount: total number of nodes
% cellCount: total number of elements
% LOC(n,1:2): localtion coordinate of node "n" in x(1) and y(2) direction
% LOCC(e,1:2): localtion coordinate of element centroid "e" in x(1) and y(2) direction

% Generate the Gauss point
GaussCount = 4;
x_gauss = zeros(cellCount, GaussCount);
y_gauss = zeros(cellCount, GaussCount);

xc_gauss = zeros(cellCount, 9*GaussCount);
yc_gauss = zeros(cellCount, 9*GaussCount);
omega_gauss = zeros(cellCount, GaussCount);

    for c=1:cellCount
    
    if GaussCount==1
        n_g = 1;
    else
        n_g = GaussCount/2;
    end
    
    x0 = LOCC(c,1)-0.5*le(1);               
    x1 = LOCC(c,1)+0.5*le(1);    
    [x_gauss_temp, x_omega_gauss] = lgwt(n_g,x0,x1);
    
    if GaussCount==4
    x_gauss(c,:) = [x_gauss_temp(1), x_gauss_temp(2), x_gauss_temp(2), x_gauss_temp(1)];   
    elseif GaussCount==1
    x_gauss(c,:) = [x_gauss_temp(1)];
    end
    
    y0 = LOCC(c,2)-0.5*le(2);               
    y1 = LOCC(c,2)+0.5*le(2);   
    [y_gauss_temp, y_omega_gauss] = lgwt(n_g,y0,y1);
    
    if GaussCount==4
    y_gauss(c,:) = [y_gauss_temp(1), y_gauss_temp(1), y_gauss_temp(2), y_gauss_temp(2)];
    elseif GaussCount==1
    y_gauss(c,:) = [y_gauss_temp(1)];
    end
    
    omega_gauss_temp = x_omega_gauss*y_omega_gauss';
    omega_gauss(c,:) = reshape(omega_gauss_temp, 1, GaussCount);
    end
    

%% Time
C                       = sqrt(E/psp);
ftime                   = 0.02;
dt                      = 0.00001;
ndt                     = round(ftime/dt) +1;
t                       = 0;

%% Boundary nodes
% Boundary coordination
x_min = 2*le(1);
x_max = 1+2*le(1);
y_min = 1*le(2);
y_max = 2*le(2);
[nfbcx,nfbcy,fbcx,fbcy]=Compute_Boundary_Nodes(nodeCount,LOC,x_max,x_min,y_max,y_min);
% nfbcx: number of boundary nodes in X direction
% nfbcy: number of boundary nodes in Y direction
% fbcx: index of all boundary nodes in X direction
% fbcy: index of all boundary nodes in Y direction

% for n=1:nodeCount
%     if LOC(n,1)<=le(1) || LOC(n,1)>=1+3*le(1)
%         nfbcy = nfbcy+1;
%         fbcy = [fbcy n];
%     end
% end

% for n=1:nodeCount
%     if LOC(n,2)<=le(2) || LOC(n,2)>=4*le(2)
%         nfbcx = nfbcx+1;
%         fbcx = [fbcx n];
%     end
% end

%% Particle generation
particle_per_cell       = 4;
spCount                 = 1*resolution*particle_per_cell;
lp(1)                   = le(1)/sqrt(particle_per_cell);                                 % size of particle in X direction
lp(2)                   = le(2)/sqrt(particle_per_cell);                                % size of particle in Y direction
x_sp                    = zeros(spCount,2);
d_sp                    = zeros(spCount,2);

sp=1;
while sp<spCount+0.0001
    for i=1:1*sqrt(particle_per_cell)
        for j=1:resolution*sqrt(particle_per_cell)
            x_sp(sp,1:2)= [2*le(2)+0.5*lp(1,1)+(j-1)*lp(1,1) 1*le(2)+0.5*lp(1,2)+(i-1)*lp(1,2)];
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
d_sp                    = zeros(spCount,2);                     % displacement
b_sp                    = zeros(spCount,2);                     % body force
s_sp                    = zeros(spCount,3);                     % Stress tensor
ds_sp                   = zeros(spCount,3);                     % Stress increment
v_ssp                   = zeros(spCount,2);                     % velocty
v_ssp1                  = zeros(spCount,2);                     % velocty
e_sp                    = zeros(spCount,3);                     % Strain tensor
de_sp                   = zeros(spCount,3);                     % Strain increment
ptraction_sp            = zeros(spCount,2);                     % traction
F_sp                    = cell(spCount,1);                      % Gradient deformation
r1_sp                   = zeros(spCount,2);
r2_sp                   = zeros(spCount,2);

b_sc                    = zeros(cellCount,2);

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

% Generate the vibration
for sp=1:spCount
    v_ssp(sp,1) = A*C*pi*sin(pi*(x_spo(sp,1)-2*le(1)));
end

% Data
sx1=[];sy1=[];sx2=[];sy2=[];
ux1=[];uy1=[];
time=[];
error_X=[];error_V=[];error_S=[];
perror_X=0;perror_V=0;perror_S=0;
dx1=[];dy1=[];
Fxx1=[];Fxx2=[];Fyy1=[];Fyy2=[];
bxx=[];byy=[];

% % start the algorithm
% % video
% timestep = 100;     % number of frame to save
% r=timestep/10;      % number of frame per second video ~200s
% 
% writerObj2           = VideoWriter('aligned1Dx.avi');
% writerObj2.FrameRate = r;    % number of frame per second
% open(writerObj2);
% 
%     for tt = 1:timestep
%     ft              = ftime/timestep*tt;
    ft=ftime;
 while t<ft+0.0001      
     t

     switch Version
         case 'MPM'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = MPM_solver(CModel,CModel_parameter,...
            nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
            nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt);
        
        case 'BSMPM'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = BSMPM_solver(CModel,CModel_parameter,...
            nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
            nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt);   

         case 'MLS_MPM'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp] = MLS_MPM_solver(CModel,CModel_parameter,...
        nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,LOCC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt);

         case 'CPDI'
        [v_ssp,x_sp,d_sp,F_sp,V_sp,s_sp,p_sp,r1_sp,r2_sp] = CPDI_solver(CModel,CModel_parameter,...
        nodeCount,spCount,cellCount,x_sp,x_spo,d_sp,le,NN,LOC,b_sp,V_sp,ptraction_sp,v_ssp,s_sp,m_sp,p_sp,...
        nfbcx,nfbcy,fbcx,fbcy,F_sp,V_spo,dt,r1_sp,r10_sp,r2_sp,r20_sp);
     end
     
 % Update time and step 
 t = t+dt;
 
     % Methods of Manufactued solution
     % Compute the analytical solutions
     for sp=1:spCount      
        % Displacement
        ux(sp) = A*sin(pi*(x_spo(sp,1)-2*le(1)))*sin(C*pi*t);        
        % Velocity
        vx(sp) = A*C*pi*sin(pi*(x_spo(sp,1)-2*le(1)))*cos(C*pi*t);        
        % Deformation gradient
        Fxx(sp) = 1 + A*pi*cos(pi*(x_spo(sp,1)-2*le(1)))*sin(C*pi*t);
        % Stress
        sx(sp) = (Lambda*log(Fxx(sp))+Mu*(Fxx(sp)^2-1))/Fxx(sp);        
        % Store analytical solutions5
        dis_ana(sp,:) =  [ux(sp)];
        vel_ana(sp,:) =  [vx(sp)];
        stress_ana(sp,:) =  [sx(sp)];
    
        % Body force applied to the particles
        b_sp(sp,1) = pi^2*ux(sp)/psp*(Lambda/Fxx(sp)^2*(1-log(Fxx(sp)))+Mu/Fxx(sp)^2*(Fxx(sp)^2+1)-E);
     end
        
    time = [time t];
    % Store for 1st particle     
    % Analytical solutions for particle 1
    ux1 = [ux1 ux(particle_index)];    sx1 = [sx1 sx(particle_index)];
    Fxx1 = [Fxx1 Fxx(particle_index)];    bxx = [bxx b_sp(particle_index,1)];
    
    % Numerical solutions for particle 1
    sx2 = [sx2 s_sp(particle_index,1)];    dx1 = [dx1 d_sp(particle_index,1)];
    Fxx2 = [Fxx2 F_sp{particle_index}(1,1)];
    
%     %  %% Plot the result
% StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,0,v_ssp,spCount,r1_sp,r2_sp);

 end

%  %% Plot the result
StressProfile1=Plot_Final(x_sp,LOC,le,d_sp,0,v_ssp,spCount,r1_sp,r2_sp);
%     frame2 = getframe(StressProfile1);
%     writeVideo(writerObj2,frame2);
%     end
%     close(writerObj2);
    
     %% Compute error
    x=0;    v=0;    s=0;
    deviation_X = zeros(spCount,1);    deviation_V = zeros(spCount,1);    deviation_S = zeros(spCount,1);
    
    % Loop all particles to store the RMS
    for sp=1:spCount      
    deviation_X(sp,1) = dis_ana(sp,1)-d_sp(sp,1);
    x=x + deviation_X(sp,:)^2;
    
    deviation_V(sp,1) = vel_ana(sp,1)-v_ssp(sp,1);
    v=v + deviation_V(sp,:)^2;
    
    deviation_S(sp,1) = stress_ana(sp,1)-s_sp(sp,1);
    s=s + deviation_S(sp,:)^2;
    end
    spCount = spCount;
    
    perror_X = sqrt(x/spCount); perror_V = sqrt(v/spCount); perror_S = sqrt(s/spCount);   
    results = [perror_X perror_V perror_S];
%     
    figure
    plot(time, ux1, time, dx1);
    title('ux');
    legend('analytical','numerical')

    figure 
    plot(time,sx1,time,sx2);
    title('stress');
    legend('analytical','numerical')
    xlabel('time')
    ylabel('stress')

    figure 
    plot(time,Fxx1,time,Fxx2);
    title('Fx');
    legend('analytical','numerical')
    
    figure
    plot(time,bxx);
    title('bx');
    legend('analytical','numerical')

    figure
    plot(time,bxx);
    title('bx');
    legend('analytical','numerical')
figure
plot(x_sp(1:resolution*3,1),stress_ana(1:resolution*3,1),x_sp(1:resolution*3,1),s_sp(1:resolution*3,1))
legend('analytical','numerical')
xlabel('position')
ylabel('stress')

figure
plot(x_sp(1:resolution*3,1),vel_ana(1:resolution*3,1),x_sp(1:resolution*3,1),v_ssp(1:resolution*3,1))
legend('analytical','numerical')
xlabel('position')
ylabel('velocity')

figure
plot(x_sp(1:resolution*3,1),dis_ana(1:resolution*3,1),x_sp(1:resolution*3,1),d_sp(1:resolution*3,1))
legend('analytical','numerical')
xlabel('position')
ylabel('displacement')