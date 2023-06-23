%% Import all relevant functions and pre-computed matrices
clear all; close all;
question1_a_function_definitions;

%% EXERCISE SOLUTION: DECENTRALIZED PROJECTED SUBGRADIENT METHOD
% Iteration parameters
max_iter = 50000;

momentum_scales = [0.7,0.8,0.9,0.95];

errors = nan(max_iter,length(momentum_scales));
error_tolerance = 1e-1;
step_size = 0.004;
% Initialize dual variables for dual subgradient iteration 
lambda_init = zeros((N-1)*n,1);
lambda = lambda_init; 
velocity = zeros(size(lambda)); % momentum variable 

% True final state
x_f_true = [-4.50061264946867 ;
            -4.07093310815196 ;
            -4.79767912346124 ;
            -4.91202943845703 ];

lb = -umax/Tfinal*ones(m*Tfinal,1);
ub = umax/Tfinal*ones(m*Tfinal,1);
options = optimoptions('quadprog','Display','off');

for momentum_index = 1:length(momentum_scales)
    momentum_scale = momentum_scales(momentum_index);
    lambda = lambda_init;
for iter=2:max_iter
    % Solve Lagrangians in distributed fashion for optimal inputs u_i    
    u_1 = quadprog(2*F_11, F_12 + lambda12(lambda)'*C_1, [],[],[],[],lb,ub,[],options);
    u_2 = quadprog(2*F_21, F_22 + (lambda23(lambda)-lambda12(lambda))'*C_2, [],[],[],[],lb,ub,[],options);
    u_3 = quadprog(2*F_31, F_32 + (lambda34(lambda)-lambda23(lambda))'*C_3, [],[],[],[],lb,ub,[],options);
    u_4 = quadprog(2*F_41, F_42 - lambda34(lambda)'*C_4, [],[],[],[],lb,ub,[],options);

    u = [u_1 ; u_2 ; u_3 ; u_4];

    % Add error to history
    errors(iter,momentum_index) = max([norm(x_f_true - x_1f(u1(u))), ...
                        norm(x_f_true - x_2f(u2(u))), ...
                        norm(x_f_true - x_3f(u3(u))), ...
                        norm(x_f_true - x_4f(u4(u)))]);

    subgradient = h(u);
    velocity = momentum_scale*velocity + step_size*subgradient;
    lambda = lambda + velocity;

    % Print messages
    fprintf("Iteration %d,  ", iter);
    fprintf("error = %g\n", errors(iter,momentum_index));

    % Convergence condition
    if abs(errors(iter,momentum_index)) < error_tolerance
        break
    end
end
end

%% PLOTTING
close all;

for index = 1:length(momentum_scales)
    semilogy(1:max_iter, errors(:,index)); hold on;
end

title("Error sequences for momentum-based accelerated subgradient method");
grid on;
xlabel("iteration"); ylabel("error");  
legend("\beta = "+momentum_scales);
ax = gca; ax.XAxis.Exponent = 0;
