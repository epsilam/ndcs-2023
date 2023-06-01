clear all;

A = [1 -1.5 ; 0 1]; B = [0 ; 1]; K = [-16/3 4];
A_K = A-B*K;

sigma = 0.5;
h_min = 0.001; % minimal sampling interval

P = sdpvar(2,2);
constraints = [P>=0, A_K'*P + P*A_K <= 0];
options = sdpsettings('verbose',0);
optimize(constraints,[],options);

Q = -value(A_K'*P + P*A_K);
P = value(P);


%% SIMULATION
xi_0 = [1;1];

steps = 1:20;

s = NaN(length(steps),1);
s(1) = 0;

state_history = NaN(length(xi_0), length(steps));
state_history(:,1) = xi_0;

for k = steps(1:end-1)
    xi_s_k = state_history(:,k);
    s(k+1) = trigger(s(k), xi_s_k, sigma,P,Q,h_min)
    state_history(:,k+1) = state_update(s(k+1), s(k), xi_s_k);
end


%% HELPER FUNCTIONS

% refer to my notes and the lecture slides for the notation. this was meant
% to keep the code as close to the mathematical interpretation as possible.

% trigger function which returns the next time s_{k+1} given s_k such that
% the performance measure specified by Q and sigma is satisfied.
function s = trigger(s_k, xi_s_k,sigma,P,Q, h_min)
    % Define nonlinear optimization constraints
    nonlcon = @(t) perf_nonlcon(t, s_k, xi_s_k,sigma,P,Q, h_min);

    s = fmincon(@(t) -t, (s_k+0.00001), [],[],[],[],[],[], nonlcon);
end

% system simulation which uses xi_s_k as the initial state (with s_k the 
% initial time) and simulates for a time length of t-s_k. returns state 
% at the end.
function xi_t = state_update(t, s_k, xi_s_k)
    A = [1 -1.5 ; 0 1]; B = [0 ; 1]; K = [-16/3 4];
    M = expm(A*(t-s_k));
    xi_t = (M - (M-eye(2))*(A\B)*K)*xi_s_k;
end

% performance measure as specified in the lecture slides.
function phi = performance_measure(xi, epsilon, sigma, P,Q)
    B = [0 ; 1]; K = [-16/3 4];

    phi = [xi' epsilon'] * [(1-sigma)*Q,  P*B*K ; (B*K)'*P, zeros(2,2)] * [xi ; epsilon];

end

% Convert the performance measure and the minimum sampling interval
% constraint into constraint inequality functions compatible with fmincon.
function [c,ceq] = perf_nonlcon(t, s_k, xi_s_k,sigma,P,Q, h_min)
    xi = @(t) state_update(t,s_k,xi_s_k);
    eps = @(t) xi_s_k - xi(t);

    c = [-performance_measure(xi(t),eps(t),sigma,P,Q); s_k-t + h_min];
    ceq = 0;
end