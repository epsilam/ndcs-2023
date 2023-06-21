%% Import all relevant functions and pre-computed matrices
clear all; close all;
question1_a_function_definitions;

%% EXERCISE SOLUTION: CENTRALIZED PROJECTED SUBGRADIENT METHOD
u = sdpvar(N*m*Tfinal,1);
Constraints = [h(u) == 0, g(u) <= 0];
Objective = [f(u)];
options = sdpsettings('verbose',0,'solver','quadprog');
optimize(Constraints,Objective,options)
u = value(u);

%% PLOTTING
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