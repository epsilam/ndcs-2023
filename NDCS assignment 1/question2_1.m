A = [1 -1.5 ; 0 1]; B = [0 ; 1];
SYS = ss(A,B,[],[]);

K1 = place(A,B,[-1+2i,-1-2i]);
K2 = place(A,B,[-1,-2]);

% Augmented static controller for delayed system, see Lecture 2 Slide 18.
K1_static = [K1 0];
K2_static = [K2 0];

dh = 0.005;
sampling_intervals = dh:dh:0.41';
n_intervals = length(sampling_intervals);

data_K1_static = NaN(1,3);
data_K2_static = NaN(1,3);
f = waitbar(0, 'Starting'); n_evals = n_intervals*(n_intervals+1)/2 - n_intervals;
for idx_h = 1:n_intervals
    h = sampling_intervals(idx_h);
    delays = sampling_intervals(1:idx_h-1);

    for idx_tau = 1:length(delays)
        tau = delays(idx_tau);

        stab = test_stability(SYS,h,tau,{K1_static,K2_static});
        data_K1_static = [data_K1_static ; h tau stab(1)];
        data_K2_static = [data_K2_static ; h tau stab(2)];
    end
    idx_waitbar = idx_h*(idx_h+1)/2 - idx_h; waitbar(idx_waitbar/n_evals, f, sprintf('Progress: %d %%', floor(idx_waitbar/n_evals*100)));
end
close(f);

figure();
subplot(2,1,1);
plot_stability(data_K1_static,"Stability of closed-loop with K1");
subplot(2,1,2);
plot_stability(data_K2_static,"Stability of closed-loop with K2");

function stability = test_stability(SYS,h,tau,controllers)
    SYS.InputDelay = tau;SYS = c2d(SYS,h);stability = NaN(length(controllers));
    for idx=1:length(controllers)
        stability(idx) = max(abs(eig(SYS.A - SYS.B*controllers{idx}))) < 1;
    end
end

function plot_stability(data_K,title_string)
    stable_h_tau_pairs  = data_K(data_K(:,3)==1,1:2);
    unstable_h_tau_pairs = data_K(data_K(:,3)==0,1:2);
    scatter(stable_h_tau_pairs(:,1),stable_h_tau_pairs(:,2),[],[0 0 0],"filled"); hold on
    scatter(unstable_h_tau_pairs(:,1),unstable_h_tau_pairs(:,2),[],[0.8 0.8 0.8],".");
    title(title_string);
    xlabel('sampling interval');
    ylabel('delay');
    xlim([0 max(data_K(:,1))]);
    ylim([0 0.2]);
    legend("stable","unstable");
end