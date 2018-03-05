function[Reconstruction] = TLS2D(n_gauss, x_gauss, y_gauss, ...
    n_particles, x_sp, le, U, x_center, Conserve)

% Reconstruct the function in each element using Taylor least square (TLS)

%% Compute cell centroids and scalling factors
x_c   = x_center(1);
y_c   = x_center(2);
x_min = x_c - 0.5*le(1);
y_min = y_c - 0.5*le(2);
x_max = x_c + 0.5*le(1);
y_max = y_c + 0.5*le(2);
dx    = le(1)/2;
dy    = le(1)/2;

%% Compute volume averages
integral_phi_4 = 1/3*(x_max^3-x_min^3)-x_c*(x_max^2-x_min^2)+...
    x_c^2*(x_max-x_min);
integral_phi_4 = integral_phi_4/(2*dx^2*(x_max-x_min));
integral_phi_5 = 1/3*(y_max^3-y_min^3)-y_c*(y_max^2-y_min^2)+...
    y_c^2*(y_max-y_min);
integral_phi_5 = integral_phi_5/(2*dy^2*(y_max-y_min));
integral_phi_6 = ((x_max^2-x_min^2)/2-x_c*(x_max-x_min))*...
    ((y_max^2-y_min^2)/2-y_c*(y_max-y_min));
integral_phi_6 = integral_phi_6/(dx*dy*(x_max-x_min)*(y_max-y_min));

%  Constract basis functions and evaluate them at the particle positions
phi_1 = ones(1, n_particles);
phi_2 = (x_sp(:,1)' - x_c*ones(1, n_particles))/dx;
phi_3 = (x_sp(:,2)' - y_c*ones(1, n_particles))/dy;
phi_4 = (x_sp(:,1)' - x_c*ones(1, n_particles)).*...
    (x_sp(:,1)' - x_c*ones(1, n_particles));
phi_4 = phi_4/(2*dx^2) - integral_phi_4*ones(1, n_particles);
phi_5 = (x_sp(:,2)' - y_c*ones(1, n_particles)).*...
    (x_sp(:,2)' - y_c*ones(1, n_particles));
phi_5 = phi_5/(2*dy^2) - integral_phi_5*ones(1, n_particles);
phi_6 = (x_sp(:,1)' - x_c*ones(1, n_particles)).*...
    (x_sp(:,2)' - y_c*ones(1, n_particles));
phi_6 = phi_6/(dx*dy) - integral_phi_6*ones(1, n_particles);

%% Constract basis functions and evaluate them at the Gauss points
phi_x_1 = ones(1, n_gauss);
phi_x_2 = (x_gauss - x_c*ones(1, n_gauss))/dx;
phi_x_3 = (y_gauss - y_c*ones(1, n_gauss))/dy;
phi_x_4 = (x_gauss - x_c*ones(1, n_gauss)).*...
    (x_gauss - x_c*ones(1, n_gauss));
phi_x_4 = phi_x_4/(2*dx^2) - integral_phi_4*ones(1, n_gauss);
phi_x_5 = (y_gauss - y_c*ones(1, n_gauss)).*...
    (y_gauss - y_c*ones(1, n_gauss));
phi_x_5 = phi_x_5/(2*dy^2) - integral_phi_5*ones(1, n_gauss);
phi_x_6 = (x_gauss - x_c*ones(1, n_gauss)).*...
    (y_gauss - y_c*ones(1, n_gauss));
phi_x_6 = phi_x_6/(dx*dy) - integral_phi_6*ones(1, n_gauss);

%% Choice of numbers of basis functions
if Conserve ~= -999 && n_particles < 3 %if we use conservation and there 
    % are less than three particles we don't need LS approximation
    Reconstruction = Conserve*ones(n_gauss,1);
else
    if Conserve == -999 
        if n_particles < 6 % use only 3 basis function
        P = [phi_1; phi_2; phi_3];
        D   = zeros (3);
        elseif n_particles >= 6 % use all 6 basis functions
        P = [phi_1; phi_2; phi_3; phi_4; phi_5; phi_6];
        D   = zeros (6);
        end
    else
        if n_particles < 6 % use only 2 basis functions
        P = [phi_2; phi_3];
        D   = zeros (2);
        elseif n_particles >= 6 % use 5 basis functions
        P = [phi_2; phi_3; phi_4; phi_5; phi_6];
        D   = zeros (5);
        end
    end

B = P;
for i = 1 : n_particles
    PP = P(:,i)*P(:,i)';
    D  = D + PP;
end

if Conserve == -999
    if n_particles < 6 %use only 3 basis function
        Phi_x = [phi_x_1', phi_x_2', phi_x_3'];
    elseif n_particles >= 6  %use all 6 basis functions
        Phi_x = [phi_x_1', phi_x_2', phi_x_3', phi_x_4', phi_x_5', phi_x_6'];
    end
    alpha = B*U;
else
    if n_particles < 6 %use only 2 basis function
        Phi_x = [phi_x_2', phi_x_3'];
    else
        Phi_x = [phi_x_2', phi_x_3', phi_x_4', phi_x_5', phi_x_6'];
    end
    alpha = B*(U-Conserve*ones(n_particles,1));
end
alpha = D\alpha;

if Conserve == -999
    Reconstruction = Phi_x*alpha;
else
    Reconstruction = Phi_x*alpha + Conserve*ones(n_gauss,1);
end

Reconstruction = Reconstruction';
end
