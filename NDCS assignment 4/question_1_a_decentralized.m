%% Import all relevant functions and pre-computed matrices
clear all; close all;
question1_a_function_definitions;

%% EXERCISE SOLUTION: DECENTRALIZED PROJECTED SUBGRADIENT METHOD
% Iteration parameters
max_iter = 50000;
errors = nan(max_iter,1);
step_size = 0.004;

% Initialize dual variables for dual subgradient iteration 
lambda_init = 1*ones((N-1)*n,1);
lambda_init = [-10.6424;-4.0134;8.7309; 2.6704; 14.7378; 15.8209; -3.0359; -6.5858; 17.9120; 17.4355; -5.5997;-5.4517];
lambda_init = [-12.5145    1.0948   10.4644   -0.8572   35.9010   42.5272  -11.7139  -19.8729   38.5687   40.5854  -11.9388  -12.5561]';
lambda_init = [-9.0624    6.4459    7.9524   -4.7764   44.9703   54.7992  -16.7634  -26.9183   45.7579   50.0472  -14.1451  -15.4598]';
%lambda_init = [-8.0542    7.8649    7.2122   -5.8195   47.1508   57.7843  -18.0332  -28.6688   47.4167   52.3028  -14.6542  -16.1520]';
lambda = lambda_init; 
u_prev = nan(N*m*Tfinal,1); % Initialize variable to store 1 previous iteration
lambda_prev = nan((N-1)*n,1);

% Keep track of final states
final_state_history_1 = nan(n,max_iter);
final_state_history_2 = nan(n,max_iter);
final_state_history_3 = nan(n,max_iter);
final_state_history_4 = nan(n,max_iter);

final_state_diffs_history = nan(N-1,max_iter);
final_state_total_diffs_history = nan(1,max_iter);

lagrangian_history = nan(max_iter,1);
%objective_function_history = nan(max_iter,1);

for iter=2:max_iter
    tic;
    % Solve Lagrangians in distributed fashion for optimal inputs u_i
    u_1 = sdpvar(m*Tfinal,1);
    u_2 = sdpvar(m*Tfinal,1);
    u_3 = sdpvar(m*Tfinal,1);
    u_4 = sdpvar(m*Tfinal,1);
    options = sdpsettings('verbose',0,'solver','quadprog');

    Obj1 = L_1(u_1,lambda); % Define objective functions for YALMIP to 
    Obj2 = L_2(u_2,lambda); % minimize, and to save the function value 
    Obj3 = L_3(u_3,lambda); % after minimization
    Obj4 = L_4(u_4,lambda);

    optimize(-umax/Tfinal <= u_1 <= umax/Tfinal, Obj1, options);
    optimize(-umax/Tfinal <= u_2 <= umax/Tfinal, Obj2, options);
    optimize(-umax/Tfinal <= u_3 <= umax/Tfinal, Obj3, options);
    optimize(-umax/Tfinal <= u_4 <= umax/Tfinal, Obj4, options);
    
    u = [value(u_1) ; value(u_2) ; value(u_3) ; value(u_4)];

    % Add newly computed quantities to history trackers
    %lagrangian_history(iter) = L(u,lambda);
    lagrangian_history(iter) = value(Obj1)+value(Obj2)+value(Obj3)+value(Obj4);
    subgradient = h(u);

    % Add final states to history
    %final_state_history_1(:,iter) = x_1f(u1(u));
    %final_state_history_2(:,iter) = x_2f(u2(u));
    %final_state_history_3(:,iter) = x_3f(u3(u));
    %final_state_history_4(:,iter) = x_4f(u4(u));

    % Add distances between final states to history
    final_state_diffs = subgradient;
    final_state_diffs_history(1,iter) = norm(final_state_diffs(0*n+1:1*n));
    final_state_diffs_history(2,iter) = norm(final_state_diffs(1*n+1:2*n));
    final_state_diffs_history(3,iter) = norm(final_state_diffs(2*n+1:3*n));

    % Update variable step size
    %if ~isnan(u_prev)
        %step_size = 0.1/norm(h(u));
        %step_size = 0.1/(1+iter);
    %end

    % Update dual variables
    lambda = lambda + step_size*subgradient; 

    sum_diffs = sum(final_state_diffs_history(:,iter));

    % Print messages
    fprintf("Iteration %d,  ", iter);
    fprintf("L - L_prev = %g,  ", lagrangian_history(iter)-lagrangian_history(iter-1));
    fprintf("L = %g,   ", lagrangian_history(iter));
    fprintf("Total x_f differences = %g\n", sum_diffs);

    toc;

    % Convergence condition
    if sum_diffs < 0.5
        break
    end
end

%% PLOTTING
figure(); 
subplot(2,1,1); hold on;
title("Dual function at each iteration of subgradient iteration");
semilogy(1:max_iter, lagrangian_history);
xlabel("iteration"); ylabel("Dual function value");
xlim([1, iter]);

subplot(2,1,2); hold on;
title("Norms of differences between final states during optimization");
semilogy(1:max_iter, final_state_diffs_history(1,:));
semilogy(1:max_iter, final_state_diffs_history(2,:));
semilogy(1:max_iter, final_state_diffs_history(3,:));
legend("x_{1f} - x_{2f}", "x_{2f} - x_{3f}", "x_{3f} - x_{4f}");
xlabel("iteration"); ylabel("norm");
xlim([1, iter]);

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