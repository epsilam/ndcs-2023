clear all;
close all;

gamma = 0.1; % stability tolerance parameter
dt = 0.01;
h_values = dt:dt:0.5;

h_tau_pairs = NaN(length(h_values)^2,2);
stability_results = NaN(length(h_values)^2,1);

index = 0;
for h=h_values
    for tau_max=dt:dt:h
        fprintf("h = %g, ", h); fprintf("tau_max = %g, ", tau_max);
        index = index + 1;
        h_tau_pairs(index,:) = [h,tau_max];
        stability_results(index) = test_stability(h,tau_max,gamma);
        fprintf("stability = %d\n", stability_results(index));
    end
end

stable_h_tau_pairs   = h_tau_pairs(stability_results==1,:);
unstable_h_tau_pairs = h_tau_pairs(stability_results==0,:);

figure(); hold on;
scatter(  stable_h_tau_pairs(:,1),  stable_h_tau_pairs(:,2),"filled","black");
scatter(unstable_h_tau_pairs(:,1),unstable_h_tau_pairs(:,2),"black");
legend("stable","LMIs found infeasible");
xlabel("h"); ylabel("\tau_{max}"); title("Stability region over (h,\tau_{max})");


function stab = test_stability(h,tau_max,gamma)
    [M0, M1, M2] = basis_elements(h, tau_max);
    [V1, V2, V3] = outer_polytopic_approximation_3(h,tau_max);

    % Closed-loop matrix outer approximations
    M_V1 = M0 + M1*V1(1) + M2*V1(2);
    M_V2 = M0 + M1*V2(1) + M2*V2(2);
    M_V3 = M0 + M1*V3(1) + M2*V3(2);
    
    P = sdpvar(3,3);
    constraints = [P >= eye(3), ...
                   M_V1' * P * M_V1 - P <= -gamma * P, ...
                   M_V2' * P * M_V2 - P <= -gamma * P, ...
                   M_V3' * P * M_V3 - P <= -gamma * P];

    options = sdpsettings('verbose',1);
    diagnostics = optimize(constraints,P,options);
    stab = strcmp('Successfully solved (LMILAB)',diagnostics.info);
end


function [V1, V2, V3] = outer_polytopic_approximation_3(h,tau_max)
    alpha_1 = @(tau) exp(h-tau);
    alpha_2 = @(tau) exp(h-tau)*(h-tau);

    V1 = [alpha_1(tau_max),alpha_2(tau_max)];
    V2 = [alpha_1(0),alpha_2(0)];
    V3 = [alpha_1(0),alpha_2(tau_max)];
end

function [M0, M1, M2] = basis_elements(h_arg, tau_arg)
    %% SYSTEM MATRICES
    syms h tau
    A = [1 -1.5 ; 0 1]; B = [0 ; 1]; K = [-16/3 4 0];
    
    F_x = expm(A*h);  F_u = (expm(A*h) - expm(A*(h-tau)))*(A\B);
    G_1 = (expm(A*(h-tau)) - eye(2))*(A\B);
    F = [F_x F_u ; zeros(1,3)]; G = [G_1 ; 1]; M = F - G*K;
    
    %% DECOMPOSE SYSTEM MATRIX INTO JORDAN BASIS
    syms alpha1 alpha2
    M = subs(M, exp(h-tau)*(h-tau), alpha2);
    M = subs(M, exp(h-tau), alpha1);
    
    % Basis elements
    M2 = diff(M, alpha2);
    M1 = diff(M, alpha1);
    M0 = M - M1*alpha1 - M2*alpha2;

    % Evaluate symbolic expressions
    h = h_arg; tau = tau_arg;
    M0 = eval(M0);
    M1 = eval(M1);
    M2 = eval(M2);
end