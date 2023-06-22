%% Import all relevant functions and pre-computed matrices
clear all; close all;
question1_a_function_definitions;

W = [0.75 0.25 0 0 ; 0.25 0.5 0.25 0 ; 0 0.25 0.5 0.25 ; 0 0 0.25 0.75];


%% EXERCISE SOLUTION: COMBINED CONCENSUS/INCREMENTAL SUBGRADIENT METHOD
% Iteration parameters
max_iter = 50000;
errors = nan(max_iter,1);
step_size = 0.002;
error_tolerance = 1e-2;

% True final state
x_f_true = [-4.50061264946867 ;
            -4.07093310815196 ;
            -4.79767912346124 ;
            -4.91202943845703 ];

% Inequality constraint bounds
lb = -umax/Tfinal*ones(m*Tfinal,1);
ub =  umax/Tfinal*ones(m*Tfinal,1);
options = optimoptions('quadprog','Display','off');

% Precompute terms for the iteration
v1 = A1^Tfinal*x01;
v2 = A2^Tfinal*x02;
v3 = A3^Tfinal*x03;
v4 = A4^Tfinal*x04;

phi = 100; % number of consensus steps;
W_phi = W^phi;
x_f_local = zeros(n,N);

for iter=2:max_iter
    tic;    
    [u_1,~,~,~,nu_1] = quadprog(2*F_11, F_12, [],[],C_1,x_f_local(:,1) - v1,lb,ub,[],options);
    [u_2,~,~,~,nu_2] = quadprog(2*F_21, F_22, [],[],C_2,x_f_local(:,2) - v2,lb,ub,[],options);
    [u_3,~,~,~,nu_3] = quadprog(2*F_31, F_32, [],[],C_3,x_f_local(:,3) - v3,lb,ub,[],options);
    [u_4,~,~,~,nu_4] = quadprog(2*F_41, F_42, [],[],C_4,x_f_local(:,4) - v4,lb,ub,[],options);

    u = [u_1 ; u_2 ; u_3 ; u_4];
    % Internal dual variables for equality constraints from quadprog
    nu = [nu_1.eqlin, nu_2.eqlin, nu_3.eqlin, nu_4.eqlin];
    
    for i=1:4
        x_if_local_terms = nan(n,N);
        for j=1:4
            x_if_local_terms(:,j) = W_phi(i,j) * (x_f_local(:,j) + step_size * nu(:,j));
        end
        x_f_local(:,i) = sum(x_if_local_terms, 2);
    end
    % Error is max distance of each agent's final state to true final state
    errors(iter) = max([norm(x_f_true - x_f_local(:,1)), ...
                  norm(x_f_true - x_f_local(:,2)), ...
                  norm(x_f_true - x_f_local(:,3)), ...
                  norm(x_f_true - x_f_local(:,4))]);

    % Print messages
    fprintf("Iteration %d,  ", iter);
    fprintf("error = %g,   ", errors(iter));

    toc;
    % Convergence condition
    if abs(errors(iter)) < error_tolerance
        break
    end
end

%% PLOTTING
close all;

figure(); 
semilogy(1:max_iter, errors);
title("Error sequence");
grid on;
xlabel("iteration"); ylabel("error");
legend("max_i |x_{i,f}^k - x_f^*|");
xlim([1, iter]); 

% subplot(3,1,3); hold on; grid on;
% title("Distance betwen true optimized final state and individual agent final state");
% plot(1:max_iter, final_state_diffs_history(:,1));
% plot(1:max_iter, final_state_diffs_history(:,2));
% plot(1:max_iter, final_state_diffs_history(:,3));
% plot(1:max_iter, final_state_diffs_history(:,4));
% legend("|x_{f} - x_{1f}|", "|x_{f} - x_{2f}|", "|x_{f} - x_{3f}|", "|x_{f} - x_{4f}|");
% xlabel("iteration"); ylabel("norm");
% xlim([1, iter]); ylim([0, 20]);

% figure();
% linear_lagrangian_step_ratios = ...
%     (f_optimal - lagrangian_history(2:end)) ./ (f_optimal - lagrangian_history(1:end-1));
% plot(1:max_iter-1, linear_lagrangian_step_ratios);
% title("Demonstration of linear convergence");
% grid on;
% xlabel("iteration"); ylabel("ratio of successive errors");
% legend("|f* - d(λ^{k+1})| / |f* - d(λ^{k})|",'FontSize',14);
% xlim([1,iter]);
