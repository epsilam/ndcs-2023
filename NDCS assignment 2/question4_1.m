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


figure(); hold on;
h=1;
alpha_1 = @(tau) exp(h-tau);
alpha_2 = @(tau) exp(h-tau)*(h-tau);
tau_values = 0:0.01:h;
plot(arrayfun(alpha_1,tau_values),arrayfun(alpha_2,tau_values),'LineWidth',2);
xlabel("\alpha_1(\tau)"); ylabel("\alpha_2(\tau)"); grid on;

V_1 = [alpha_1(tau_values(end)),alpha_2(tau_values(end))];
V_2 = [alpha_1(tau_values(1)),alpha_2(tau_values(1))];
V_3 = [alpha_1(tau_values(1)),alpha_2(tau_values(end))];
V = [V_1 ; V_2 ; V_3];

K = convhull(V(:,1),V(:,2));
plot(V(K,1),V(K,2),'LineWidth',2)
title("Polytopic outer approximation of (\alpha_1(\tau),\alpha_2(\tau)) trajectory with h=1");
legend("(\alpha_1(\tau),\alpha_2(\tau)) trajectory through \tau\in [0,h])", ...
       "polytopic outer approximation")