A = [1 -1.5 ; 0 1]; B = [0 ; 1];
SYS = ss(A,B,[],[]);
K1 = place(A,B,[-1+2i,-1-2i]);

h = 0.25; tau=0.08;
F_x = expm(A*h);
            F_u = (expm(A*h) - expm(A*(h-tau)))*(A\B);
            G_1 = (expm(A*(h-tau)) - eye(2))*(A\B);
            A_aug = [F_x        F_u zeros(2,1); 
                     zeros(1,2) 0   0         ;
                     zeros(1,2) 1   0         ];
            B_aug = [G_1 ; 1 ; 0];
K1_dynamic = place(A_aug,B_aug,exp([-1+2i,-1-2i,-2,-5]));


% Augmented static controller for delayed system, see Lecture 2 Slide 18.
K1_static = [K1 0];

% Define range of sampling intervals h>0 to test for closed-loop stability
dt                 = 0.005;
hmax               = 0.4;
sampling_intervals = 0.01:dt:hmax;
delay_times        = 0.01:dt:1.5*hmax;

f = waitbar(0, 'Starting');
n_intervals = length(sampling_intervals);
n_evals = n_intervals*(n_intervals+1)/2 - n_intervals;

data_K1_dynamic = NaN(1,3);

for idx_h = 1:n_intervals
    h = sampling_intervals(idx_h);
    delays = 0.01:dt:1.5*h;
    for idx_tau = 1:(length(delays)-1)
        tau = delays(idx_tau);
        stab = test_stability(SYS,h,tau,{K1_dynamic});
        data_K1_dynamic = [data_K1_dynamic ; h tau stab(1)];
    end

    idx_waitbar = idx_h*(idx_h+1)/2 - idx_h;
    waitbar(idx_waitbar/n_evals, f, sprintf('Progress: %d %%', floor(idx_waitbar/n_evals*100)));
end
close(f);

[stab_K1_dynamic, unstab_K1_dynamic] = stable_unstable_pairs(data_K1_dynamic);

%%
figure;
scatter(stab_K1_dynamic(:,1),stab_K1_dynamic(:,2),[],"black","o","filled"); hold on;
scatter(unstab_K1_dynamic(:,1),unstab_K1_dynamic(:,2),1,[0.8,0.8,0.8],"."); hold on;
plot(0:0.1:1,0:0.1:1,"red",'LineWidth',2);
title("Stability of closed-loop system with controller K1");
xlabel('sampling interval');ylabel('delay');
    xlim([0 0.4]);
    ylim([0 0.2]);
legend("stable","unstable","h=tau line");

function stability = test_stability(SYS,h,tau,controllers)
    %SYS.InputDelay = tau;SYS = c2d(SYS,h);stability = NaN(length(controllers));
    A = SYS.A;
    B = SYS.B;
    F_x = expm(A*h);
    for idx=1:length(controllers)
        if tau <= h
            F_u = (expm(A*h) - expm(A*(h-tau)))*(A\B);
            G_1 = (expm(A*(h-tau)) - eye(2))*(A\B);
            A_aug = [F_x        F_u zeros(2,1); 
                     zeros(1,2) 0   0         ;
                     zeros(1,2) 1   0         ];
            B_aug = [G_1 ; 1 ; 0];
            K = controllers{idx};
            stability(idx) = max(abs(eig(A_aug - B_aug*K))) < 1;
        elseif tau>h
            F_u1 = (expm(A*(2*h-tau)) - eye(2))*(A\B);
            F_u2 = (expm(A*h) - expm(A*(2*h-tau)))*(A\B);
            A_aug = [F_x        F_u1 F_u2 ; 
                     zeros(1,2) 0    0    ;
                     zeros(1,2) 1    0    ];
            B_aug = [zeros(2,1) ; 1 ; 0];
            K = controllers{idx};
            stability(idx) = max(abs(eig(A_aug - B_aug*K))) < 1;
        end

    end
end

function [stable_h_tau_pairs, unstable_h_tau_pairs] = stable_unstable_pairs(data_K)
    stable_h_tau_pairs  = data_K(data_K(:,3)==1,1:2);
    unstable_h_tau_pairs = data_K(data_K(:,3)==0,1:2);
end