% Note to reader: all notation is designed to match the mathematical
% notation given in the report as closely as possible. 

%% Initialize system matrices and standard variables
aircraft;

%% Functions to extract distributed components from large vectors
% Input vector
ui_dim = m*Tfinal;
u1 = @(u) u(0*ui_dim+1 : 1*ui_dim);
u2 = @(u) u(1*ui_dim+1 : 2*ui_dim);
u3 = @(u) u(2*ui_dim+1 : 3*ui_dim);
u4 = @(u) u(3*ui_dim+1 : 4*ui_dim);

% Dual variables for equality constraints
lambda12 = @(lambda) lambda(0*n+1 : 1*n);
lambda23 = @(lambda) lambda(1*n+1 : 2*n);
lambda34 = @(lambda) lambda(2*n+1 : 3*n);

%% Define minimization objective, constraint functions, and duals
% Precompute for f (minimization objective function)
[F_11, F_12, F_13] = f_matrices(1); f_1 = @(u_1) u_1'*F_11*u_1 + F_12*u_1 + F_13;
[F_21, F_22, F_23] = f_matrices(2); f_2 = @(u_2) u_2'*F_21*u_2 + F_22*u_2 + F_23;
[F_31, F_32, F_33] = f_matrices(3); f_3 = @(u_3) u_3'*F_31*u_3 + F_32*u_3 + F_33;
[F_41, F_42, F_43] = f_matrices(4); f_4 = @(u_4) u_4'*F_41*u_4 + F_42*u_4 + F_43;
f = @(u) f_1(u1(u)) + f_2(u2(u)) + f_3(u3(u)) + f_4(u4(u));

% Precompute for h (equality constraint function)
zero = zeros(n,m*Tfinal);
A_h = [C(1,Tfinal), -C(2,Tfinal),         zero,        zero ; 
              zero,  C(2,Tfinal), -C(3,Tfinal),        zero ;
              zero,         zero,  C(3,Tfinal), -C(4,Tfinal)];
b_h = [A1^Tfinal*x01 - A2^Tfinal*x02 ; 
       A2^Tfinal*x02 - A3^Tfinal*x03 ; 
       A3^Tfinal*x03 - A4^Tfinal*x04];
h = @(u) A_h * u + b_h;

% Define g function (inequality constraint function)
g = @(u) [g_i(u1(u)) ; g_i(u2(u)) ; g_i(u3(u)) ; g_i(u4(u))];

% Define final state functions
x_1f = @(u_1) A1^Tfinal*x01 + C(1,Tfinal)*u_1;
x_2f = @(u_2) A2^Tfinal*x02 + C(2,Tfinal)*u_2;
x_3f = @(u_3) A3^Tfinal*x03 + C(3,Tfinal)*u_3;
x_4f = @(u_4) A4^Tfinal*x04 + C(4,Tfinal)*u_4;

% Define distributed Lagrangians
L_1 = @(u_1,lambda) f_1(u_1)                      + lambda12(lambda)'*x_1f(u_1);
L_2 = @(u_2,lambda) f_2(u_2) + (lambda23(lambda) - lambda12(lambda))'*x_2f(u_2);
L_3 = @(u_3,lambda) f_3(u_3) + (lambda34(lambda) - lambda23(lambda))'*x_3f(u_3);
L_4 = @(u_4,lambda) f_4(u_4)                      - lambda34(lambda)'*x_4f(u_4);
L = @(u,lambda) L_1(u1(u),lambda) ...
                 + L_2(u2(u),lambda) ...
                 + L_3(u3(u),lambda) ...
                 + L_4(u4(u),lambda);



%% HELPER FUNCTIONS
function [A, B, x0, n, m, N, umax, Tfinal] = select_system(i)
    aircraft;
    if     i==1; A=A1; B=B1; x0=x01;
    elseif i==2; A=A2; B=B2; x0=x02;
    elseif i==3; A=A3; B=B3; x0=x03;
    elseif i==4; A=A4; B=B4; x0=x04;
    end
end

function [F_i1, F_i2, F_i3] = f_matrices(i)
    % Creates matrices F_i1, F_i2, F_i3 so that f_i(u) is the quadratic
    % f_i(u) = u'*F_i1*u + F_i2*u + F_i3
    [A, B, x0, n, m, N, umax, Tfinal] = select_system(i);

    F_i1 = eye(m*Tfinal);
    F_i2 = zeros(n,m*Tfinal);
    F_i3 = zeros(n,n);
    for t=1:(Tfinal-1)
        F_i1 = F_i1 + C(i,t)'*C(i,t);
        F_i2 = F_i2 + (A^t)'*C(i,t);
        F_i3 = F_i3 + (A')^t*A^t;
    end
    F_i2 = 2*x0'*F_i2;
    F_i3 = x0'*F_i3*x0;
end

% Reversed controllability matrix for precomputing state trajectories
function out = C(i,t)
    [A, B, x0, n, m, N, umax, Tfinal] = select_system(i);

    out = zeros(n,m*Tfinal);
    for t_index=0:t-1
        out(:, m*t_index+1:m*(t_index+1)) = A^(t-1-t_index)*B;
    end
end

% Distributed g function
function out = g_i(u_i)
    Tfinal=5;
    umax=60;
    out = [u_i ; -u_i] - umax/Tfinal;
end

function plot_system_trajectories(u)
    aircraft;
    ui_dim = m*Tfinal;
    u_1 = u(0*ui_dim+1 : 1*ui_dim);
    u_2 = u(1*ui_dim+1 : 2*ui_dim);
    u_3 = u(2*ui_dim+1 : 3*ui_dim);
    u_4 = u(3*ui_dim+1 : 4*ui_dim);

    x_1 = sim_system(1,u_1);
    x_2 = sim_system(2,u_2);
    x_3 = sim_system(3,u_3);
    x_4 = sim_system(4,u_4);

    figure();
    % Component 1 of each system (aircraft x position)
    subplot(2,2,1); hold on;
        plot(0:Tfinal,x_1(1,:));
        plot(0:Tfinal,x_2(1,:));
        plot(0:Tfinal,x_3(1,:));
        plot(0:Tfinal,x_4(1,:));
        legend("Aircraft 1","Aircraft 2","Aircraft 3","Aircraft 4");
        title("Aircraft trajectory x position");
    
    % Component 2 of each system (aircraft y position)
    subplot(2,2,1);
        plot(0:Tfinal,x_1(2,:));
        plot(0:Tfinal,x_2(2,:));
        plot(0:Tfinal,x_3(2,:));
        plot(0:Tfinal,x_4(2,:));
        legend("Aircraft 1","Aircraft 2","Aircraft 3","Aircraft 4");
        title("Aircraft trajectory y position");

    % Component 3 of each system (aircraft dx/dt position)
    subplot(2,2,1);
        plot(0:Tfinal,x_1(3,:));
        plot(0:Tfinal,x_2(3,:));
        plot(0:Tfinal,x_3(3,:));
        plot(0:Tfinal,x_4(3,:));
        legend("Aircraft 1","Aircraft 2","Aircraft 3","Aircraft 4");
        title("Aircraft trajectory x velocity");

    % Component 4 of each system (aircraft dy/dt position)
    subplot(2,2,1);
        plot(0:Tfinal,x_1(1,:));
        plot(0:Tfinal,x_2(1,:));
        plot(0:Tfinal,x_3(1,:));
        plot(0:Tfinal,x_4(1,:));
        legend("Aircraft 1","Aircraft 2","Aircraft 3","Aircraft 4");
        title("Aircraft trajectory y velocity");
end

function x = sim_system(i,u_i)
    [A, B, x0, n, m, N, umax, Tfinal] = select_system(i);
    x = nan(n,Tfinal+1);
    x(:,1)=x0;
    u_reshaped = reshape(u_i,[m,Tfinal]);
    for t=1:Tfinal
        x(:,t+1) = A*x(:,t) + B*u_reshaped(:,t);
    end
end