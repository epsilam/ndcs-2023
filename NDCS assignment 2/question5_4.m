close all; clear all;

h_avg = 0.4675;
h_PETC = h_avg/2;

A = [1 -1.5 ; 0 1]; B = [0 ; 1]; K = [-16/3 4];
A_K = A-B*K;

P = sdpvar(2,2);
constraints = [P>=0, A_K'*P + P*A_K <= 0];
options = sdpsettings('verbose',0);
optimize(constraints,[],options);

Q = -value(A_K'*P + P*A_K);
P = value(P);

sigma_values = [0.0001, 0.001, 0.01, 0.2, 0.6, 0.9, 0.99];
initial_conditions = [0.2 0.15; 1 2 ; 4 5 ; 7 2]';
steps = 1:1000;
t_upperlimit = 5;

% for each value of sigma and each initial condition, we save how many
% timesteps were needed to stay within the performance bound on the time
% interval [0,1] seconds
num_timesteps_computed = NaN(length(sigma_values),size(initial_conditions,2));

avg_intersample_times  = NaN(length(sigma_values),size(initial_conditions,2));

colors = [1 0 0 ; 
          0 1 0 ; 
          0 0 1 ; 
          0.5 0.5 0 ; 
          0.5 0 0.5 ;
          0 0.5 0.5 ; 
          0.5 0.5 0.5];

figure(); hold on;
for sigma_index=1:length(sigma_values)

    sigma = sigma_values(sigma_index)
    color = colors(sigma_index,:);


    for sim_index = 1:size(initial_conditions,2)
        xi_0 = initial_conditions(:,sim_index);
        
        % Array of varying sample times s_k
        s = NaN(length(steps),1);
        s(1) = 0;
        
        state_history = NaN(length(xi_0), length(steps));
        state_history(:,1) = xi_0;
        
        % Simulate with event-triggered control
        for k = steps(1:end-1)
            xi_s_k = state_history(:,k);
            s_k1 = trigger_PETC(s(k), xi_s_k, sigma,P,Q,h_PETC);
            if s_k1 > t_upperlimit
                break
            else
                s(k+1) = s_k1;
                state_history(:,k+1) = state_update(s(k+1), s(k), xi_s_k);
            end
        end
        num_timesteps_computed(sigma_index,sim_index) = length(s(~isnan(s)));
        
        % Plot state trajectories
        if sim_index == 1
            plot(s,vecnorm(state_history,2),'-o', 'Color',color, 'MarkerFaceColor',color, "DisplayName",sprintf("sigma = %g", sigma));
        else
            plot(s,vecnorm(state_history,2),'-o', 'Color',color, 'MarkerFaceColor',color, "DisplayName",'');
        end
        
        s = s(~isnan(s));
        intersample_times = s(2:end) - s(1:end-1);
        avg_intersample_times(sigma_index,sim_index) = mean(intersample_times);
    end
end
xlabel("t"); ylabel("|\xi(t)|");
title("Trajectories for various initial conditions and sigma with PETC");
legend();
mean(num_timesteps_computed,2)
mean(avg_intersample_times,2)
hold off;

figure();
plot(sigma_values,mean(num_timesteps_computed,2),'-o'); xlabel("\sigma"); ylabel("timesteps");
title("Average number of timesteps in interval [0,1] seconds for various \sigma");




function s = trigger_PETC(s_k, xi_s_k,sigma,P,Q, h_PETC)
    % search through multiples of h_PETC
    
    for r=1:1000
        t = s_k + r*h_PETC;
        xi = state_update(t, s_k, xi_s_k);
        eps = xi_s_k - xi;
        if performance_measure(xi, eps, sigma, P,Q) < 0
            break
        end
    end
    s = t;
end

function xi_t = state_update(t, s_k, xi_s_k)
    A = [1 -1.5 ; 0 1]; B = [0 ; 1]; K = [-16/3 4];
    M = expm(A*(t-s_k));
    xi_t = (M - (M-eye(2))*(A\B)*K)*xi_s_k;
end

function phi = performance_measure(xi, epsilon, sigma, P,Q)
    B = [0 ; 1]; K = [-16/3 4];
    phi = [xi' epsilon'] * [(1-sigma)*Q,  P*B*K ; (B*K)'*P, zeros(2,2)] * [xi ; epsilon];
end