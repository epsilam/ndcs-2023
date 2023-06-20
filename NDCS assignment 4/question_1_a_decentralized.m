%% Import all relevant functions and pre-computed matrices
clear all; close all;
question1_a_function_definitions;

%% EXERCISE SOLUTION: DECENTRALIZED PROJECTED SUBGRADIENT METHOD
% Iteration parameters
max_iter = 50000;
errors = nan(max_iter,1);
step_size = 0.004;

% Initialize dual variables for dual subgradient iteration 
lambda_init = zeros((N-1)*n,1);
lambda = lambda_init; 

% Keep track of final states
final_state_history_1 = nan(n,max_iter);
final_state_history_2 = nan(n,max_iter);
final_state_history_3 = nan(n,max_iter);
final_state_history_4 = nan(n,max_iter);

final_state_diffs_history = nan(max_iter,N);
final_state_total_diffs_history = nan(max_iter,1);

lagrangian_history = nan(max_iter,1);

% True final state
x_f_true = [-4.5006 ; -4.0709 ; -4.7977 ; -4.9120];

% optimal objective function value
f_optimal = 2252.44178420712;

lb = -umax/Tfinal*ones(m*Tfinal,1);
ub = umax/Tfinal*ones(m*Tfinal,1);
options = optimoptions('quadprog','Display','off');

for iter=2:max_iter
    tic;
    % Solve Lagrangians in distributed fashion for optimal inputs u_i    
    u_1 = quadprog(2*F_11, F_12 + lambda12(lambda)'*C_1, [],[],[],[],lb,ub,[],options);
    u_2 = quadprog(2*F_21, F_22 + (lambda23(lambda)-lambda12(lambda))'*C_2, [],[],[],[],lb,ub,[],options);
    u_3 = quadprog(2*F_31, F_32 + (lambda34(lambda)-lambda23(lambda))'*C_3, [],[],[],[],lb,ub,[],options);
    u_4 = quadprog(2*F_41, F_42 - lambda34(lambda)'*C_4, [],[],[],[],lb,ub,[],options);

    u = [u_1 ; u_2 ; u_3 ; u_4];

    % Add newly computed quantities to history trackers
    lagrangian_history(iter) = L(u,lambda);
    subgradient = h(u);

    % Add distances between final states to history
    final_state_diffs_history(iter,:) = [norm(x_f_true - x_1f(u1(u))), ...
                                         norm(x_f_true - x_2f(u2(u))), ...
                                         norm(x_f_true - x_3f(u3(u))), ...
                                         norm(x_f_true - x_4f(u4(u)))];

    % Update variable step size
    %if true
        %step_size = 0.1/norm(h(u));
        %step_size = 1/(5+iter);
    %end

    % Update dual variables
    lambda = lambda + step_size*subgradient; 

    % Define print quantities
    error = f_optimal - lagrangian_history(iter);
    %sum_diffs = sum(final_state_diffs_history(iter,:));
    

    % Print messages
    fprintf("Iteration %d,  ", iter);
    fprintf("f* - L = %g,   ", error);
    %fprintf("Total x_f differences = %g\n", sum_diffs);

    toc;

    % Convergence condition
    if error < 1e-3
        break
    end
end

%% PLOTTING
close all;

figure(); 
subplot(3,1,1); 
plot(1:max_iter, lagrangian_history);
grid on;
title("Dual function at each iteration of subgradient iteration");
xlabel("iteration"); ylabel("Dual function value");
xlim([1, iter]);

subplot(3,1,2); 
semilogy(1:max_iter, f_optimal - lagrangian_history);
title("Error sequence (diff. between optimal objective value and dual)");
grid on;
xlabel("iteration"); ylabel("Dual objective error");
legend("f^* - d(\lambda)");
xlim([1, iter]); 

subplot(3,1,3); hold on; grid on;
title("Distance betwen true optimized final state and individual agent final state");
plot(1:max_iter, final_state_diffs_history(:,1));
plot(1:max_iter, final_state_diffs_history(:,2));
plot(1:max_iter, final_state_diffs_history(:,3));
plot(1:max_iter, final_state_diffs_history(:,4));
legend("|x_{f} - x_{1f}|", "|x_{f} - x_{2f}|", "|x_{f} - x_{3f}|", "|x_{f} - x_{4f}|");
xlabel("iteration"); ylabel("norm");
xlim([1, iter]); ylim([0, 20]);

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