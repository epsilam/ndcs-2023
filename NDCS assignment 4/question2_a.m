%% Import all relevant functions and pre-computed matrices
clear all; close all;
question1_a_function_definitions;

%% EXERCISE SOLUTION: ADMM
% Iteration parameters
max_iter = 5000;

rho_list = [0.01,0.02,0.05,0.1,0.2,0.5,1,5];

errors = nan(max_iter,length(rho_list));
error_tolerance = 1e-3;


% Initialize dual variables for dual subgradient iteration 
lambda_init = zeros((N-1)*n,1);
lambda = lambda_init; 

% True final state
x_f_true = [-4.50061264946867 ;
            -4.07093310815196 ;
            -4.79767912346124 ;
            -4.91202943845703 ];


lb = -umax/Tfinal*ones(m*Tfinal,1);
ub = umax/Tfinal*ones(m*Tfinal,1);
options = optimoptions('quadprog','Display','off');


for rho_index = 1:length(rho_list)
    rho = rho_list(rho_index);
    % Precompute quantities
    H_11 = 2*F_11 + rho*(C_1'*C_1);
    H_21 = 2*F_21 + rho*(C_2'*C_2);
    H_31 = 2*F_31 + rho*(C_3'*C_3);
    H_41 = 2*F_41 + rho*(C_4'*C_4);
    
    % Initialize dual variables
    lambda_1 = zeros(n,1);
    lambda_2 = zeros(n,1);
    lambda_3 = zeros(n,1);
    lambda_4 = zeros(n,1);
    
    % Initialize ADMM estimate of final state 
    x_f = zeros(n,1);
for iter=2:max_iter
    % Minimize augmented Lagrangians
    H_12 = F_12 + (lambda_1 + rho*(v_1 - x_f))'*C_1;
    H_22 = F_22 + (lambda_2 + rho*(v_2 - x_f))'*C_2;
    H_32 = F_32 + (lambda_3 + rho*(v_3 - x_f))'*C_3;
    H_42 = F_42 + (lambda_4 + rho*(v_4 - x_f))'*C_4;
    u_1 = quadprog(H_11, H_12, [],[],[],[],lb,ub,[],options);
    u_2 = quadprog(H_21, H_22, [],[],[],[],lb,ub,[],options);
    u_3 = quadprog(H_31, H_32, [],[],[],[],lb,ub,[],options);
    u_4 = quadprog(H_41, H_42, [],[],[],[],lb,ub,[],options);

    u = [u_1 ; u_2 ; u_3 ; u_4];

    x_f = 1/4 * sum([C_1*u_1 + v_1 + lambda_1/rho, ... 
                     C_2*u_2 + v_2 + lambda_2/rho, ...
                     C_3*u_3 + v_3 + lambda_3/rho, ...
                     C_4*u_4 + v_4 + lambda_4/rho], 2);

    lambda_1 = lambda_1 + rho*(C_1*u_1 + v_1 - x_f);
    lambda_2 = lambda_2 + rho*(C_2*u_2 + v_2 - x_f);
    lambda_3 = lambda_3 + rho*(C_3*u_3 + v_3 - x_f);
    lambda_4 = lambda_4 + rho*(C_4*u_4 + v_4 - x_f);


    % Add error to history
    errors(iter,rho_index) = norm(x_f - x_f_true);

    % Print messages
    fprintf("Iteration %d,  ", iter);
    fprintf("error = %g\n", errors(iter,rho_index));
    % Convergence condition
    if abs(errors(iter,rho_index)) < error_tolerance
        break
    end
end
end

%% PLOTTING
close all;

figure();
for rho_index = 1:length(rho_list)
    semilogy(1:max_iter, errors(:,rho_index)); hold on;
end

title("Error sequences");
grid on;
xlabel("iteration"); ylabel("error");
legend("\rho = "+rho_list);
ylim([1e-3, 1e1]);
xlim([1, 2000]); 






plot_system_trajectories(u);


% Auxillary plotting functions
function [A, B, x0, n, m, N, umax, Tfinal] = select_system(i)
    aircraft;
    if     i==1; A=A1; B=B1; x0=x01;
    elseif i==2; A=A2; B=B2; x0=x02;
    elseif i==3; A=A3; B=B3; x0=x03;
    elseif i==4; A=A4; B=B4; x0=x04;
    end
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
        xlabel("time"); ylabel("x");
    
    % Component 2 of each system (aircraft y position)
    subplot(2,2,2); hold on;
        plot(0:Tfinal,x_1(2,:));
        plot(0:Tfinal,x_2(2,:));
        plot(0:Tfinal,x_3(2,:));
        plot(0:Tfinal,x_4(2,:));
        legend("Aircraft 1","Aircraft 2","Aircraft 3","Aircraft 4");
        title("Aircraft trajectory y position");
        xlabel("time"); ylabel("y");

    % Component 3 of each system (aircraft dx/dt position)
    subplot(2,2,3); hold on;
        plot(0:Tfinal,x_1(3,:));
        plot(0:Tfinal,x_2(3,:));
        plot(0:Tfinal,x_3(3,:));
        plot(0:Tfinal,x_4(3,:));
        legend("Aircraft 1","Aircraft 2","Aircraft 3","Aircraft 4");
        title("Aircraft trajectory x velocity");
        xlabel("time"); ylabel("dx/dt");

    % Component 4 of each system (aircraft dy/dt position)
    subplot(2,2,4); hold on;
        plot(0:Tfinal,x_1(1,:));
        plot(0:Tfinal,x_2(1,:));
        plot(0:Tfinal,x_3(1,:));
        plot(0:Tfinal,x_4(1,:));
        legend("Aircraft 1","Aircraft 2","Aircraft 3","Aircraft 4");
        title("Aircraft trajectory y velocity");
        xlabel("time"); ylabel("dy/dt");
end