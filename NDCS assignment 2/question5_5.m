sigma = 0.2;
steps = 1:1000;
t_upperlimit = 5;
xi_0 = [1 ; 2];
h_avg = 0.4675;
h_PETC = h_avg/2;
h_min_ETC = 0.001;

A = [1 -1.5 ; 0 1]; B = [0 ; 1]; K = [-16/3 4]; A_K = A-B*K;
P = sdpvar(2,2);
constraints = [P>=0, A_K'*P + P*A_K <= 0];
options = sdpsettings('verbose',0);
optimize(constraints,[],options);
Q = -value(A_K'*P + P*A_K);
P = value(P);

figure(); hold on;
simulate_PETC(sigma,t_upperlimit,steps,xi_0,h_PETC,   P,Q);
simulate_ETC( sigma,t_upperlimit,steps,xi_0,h_min_ETC,P,Q);
legend("PETC","ETC");



function simulate_PETC(sigma,t_upperlimit,steps,xi_0,h_PETC,P,Q)
    state_history = NaN(length(xi_0), length(steps));
    state_history(:,1) = xi_0;
    s = NaN(length(steps),1);
    s(1) = 0;
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
    plot(s,vecnorm(state_history,2),'-o');
end

function simulate_ETC(sigma,t_upperlimit,steps,xi_0,h_min,P,Q)
    state_history = NaN(length(xi_0), length(steps));
    state_history(:,1) = xi_0;
    s = NaN(length(steps),1); 
    s(1) = 0;
    for k = steps(1:end-1)
        xi_s_k = state_history(:,k);
        s_k1 = trigger_ETC(s(k), xi_s_k, sigma,P,Q,h_min);
        if s_k1 > t_upperlimit
            break
        else
            s(k+1) = s_k1;
            state_history(:,k+1) = state_update(s(k+1), s(k), xi_s_k);
        end
    end
    plot(s,vecnorm(state_history,2),'-o');
end

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

function s = trigger_ETC(s_k, xi_s_k,sigma,P,Q, h_min)
    % Define nonlinear optimization constraints
    nonlcon = @(t) perf_nonlcon(t, s_k, xi_s_k,sigma,P,Q, h_min);
    
    opts = optimoptions(@fmincon, 'Display','off');
    s = fmincon(@(t) -t, (s_k+0.000001), [],[],[],[],[],[], nonlcon, opts);
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

function [c,ceq] = perf_nonlcon(t, s_k, xi_s_k,sigma,P,Q, h_min)
    xi = @(t) state_update(t,s_k,xi_s_k);
    eps = @(t) xi_s_k - xi(t);
    c = [-performance_measure(xi(t),eps(t),sigma,P,Q); s_k-t + h_min];
    ceq = 0;
end