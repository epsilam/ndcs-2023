%% Import all relevant functions and pre-computed matrices
clear all; close all;
question1_a_function_definitions;

W = [0.75 0.25 0 0 ; 0.25 0.5 0.25 0 ; 0 0.25 0.5 0.25 ; 0 0 0.25 0.75];

%% EXERCISE SOLUTION: COMBINED CONCENSUS/INCREMENTAL SUBGRADIENT METHOD
% Iteration parameters
max_iter = 7000;

phi_list = [1, 2, 5, 10, 20, 50, 100, 1000];

errors = nan(max_iter,length(phi_list));
step_size = 0.004;
error_tolerance = 1e-10;

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


for phi_index = 1:length(phi_list)
phi = phi_list(phi_index); % number of consensus steps;
W_phi = W^phi;
x_f_local = zeros(n,N);

for iter=2:max_iter  
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
    errors(iter,phi_index) = max([norm(x_f_true - x_f_local(:,1)), ...
                  norm(x_f_true - x_f_local(:,2)), ...
                  norm(x_f_true - x_f_local(:,3)), ...
                  norm(x_f_true - x_f_local(:,4))]);

    % Print messages
    fprintf("Iteration %6d,  ", iter);
    fprintf("error = %2.8f\n", errors(iter,phi_index));
    % Convergence condition
    if abs(errors(iter,phi_index)) < error_tolerance
        break
    end
end
end

%% PLOTTING
close all;

figure(); 

for phi_index = 1:length(phi_list)
    semilogy(1:max_iter, errors(:,phi_index)); hold on;
end

title("Error sequence");
grid on;
xlabel("iteration"); ylabel("error");
legend("\phi = "+phi_list);
xlim([1, 7000]); 


figure();
errors_shifted = nan(size(errors));
for phi_index = 1:length(phi_list)-2
    errors_shifted(:,phi_index) = errors(:,phi_index) - min(errors(:,phi_index));
    semilogy(1:max_iter, errors_shifted(:,phi_index)); hold on;
end

title("Error sequences shifted to converge to zero");
grid on;
xlabel("iteration"); ylabel("shifted error");
legend("\phi = "+phi_list);
xlim([1, 7000]); 
