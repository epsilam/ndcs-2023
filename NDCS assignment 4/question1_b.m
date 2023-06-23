%% Import all relevant functions and pre-computed matrices
clear all; close all;
question1_a_function_definitions;

% Iteration parameters
max_iter = 20000;
step_sizes = [0.0038, 0.004, 0.0042, 0.0044, 0.0046, 0.0047, 0.00472, 0.00474];
step_coeffs_variable = [0.5, 0.2, 0.1];
error_tolerance = 1e-1;

errors = nan(max_iter,length(step_sizes));

% For each constant step size, record the num of iterations to convergence
num_iterations_to_convergence = nan(1,length(step_sizes));

% Initialize dual variables for dual subgradient iteration 
lambda_init = zeros((N-1)*n,1);

% True final state
x_f_true = [-4.50061264946867 ;
            -4.07093310815196 ;
            -4.79767912346124 ;
            -4.91202943845703 ];

lb = -umax/Tfinal*ones(m*Tfinal,1);
ub = umax/Tfinal*ones(m*Tfinal,1);
options = optimoptions('quadprog','Display','off');

for step_size_index = 1:length(step_coeffs_variable)
    %step_size_coeff = step_sizes(step_size_index);
    %step_size = step_size_coeff;
    step_size_coeff = step_coeffs_variable(step_size_index);
    lambda = lambda_init; 

    for iter=1:max_iter
        % Solve Lagrangians in distributed fashion for optimal inputs u_i    
        u_1 = quadprog(2*F_11, F_12 + lambda12(lambda)'*C_1, [],[],[],[],lb,ub,[],options);
        u_2 = quadprog(2*F_21, F_22 + (lambda23(lambda)-lambda12(lambda))'*C_2, [],[],[],[],lb,ub,[],options);
        u_3 = quadprog(2*F_31, F_32 + (lambda34(lambda)-lambda23(lambda))'*C_3, [],[],[],[],lb,ub,[],options);
        u_4 = quadprog(2*F_41, F_42 - lambda34(lambda)'*C_4, [],[],[],[],lb,ub,[],options);
        u = [u_1 ; u_2 ; u_3 ; u_4];
        
        % Update variable step size
        if true
            step_size = 0.003*(1 + 1/iter^step_size_coeff);
        end
    
        % Update dual variables
        subgradient = h(u);
        lambda = lambda + step_size*subgradient; 
    
        % Add error to history
        errors(iter,step_size_index) = max([norm(x_f_true - x_1f(u1(u))), ...
                                            norm(x_f_true - x_2f(u2(u))), ...
                                            norm(x_f_true - x_3f(u3(u))), ...
                                            norm(x_f_true - x_4f(u4(u)))]);

        % Print messages
        fprintf("Iteration %d,  ", iter);
        fprintf("error = %2.8f\n%", errors(iter, step_size_index));
    
        % Convergence condition
        if abs(errors(iter,step_size_index)) < error_tolerance
            break
        end
    end
    num_iterations_to_convergence(step_size_index) = iter;
end

%% PLOTTING
figure();
s
% subplot(2,1,1);
% for step_size_index = 1:length(step_sizes)
%     semilogy(1:max_iter, errors(:,step_size_index)); hold on;
% end
% 
% title("Error sequences for constant step sizes");
% grid on;
% xlabel("iteration"); ylabel("error");
% legend("\alpha^k = "+step_sizes);

for step_size_index = 1:length(step_coeffs_variable)
    semilogy(1:max_iter, errors(:,step_size_index)); hold on;
end

title("Error sequences for variable step sizes");
grid on;
xlabel("iteration"); ylabel("error");  
legend("\alpha^k = 0.003(1+1/k^"+step_coeffs_variable+")");
ax = gca; ax.XAxis.Exponent = 0;
