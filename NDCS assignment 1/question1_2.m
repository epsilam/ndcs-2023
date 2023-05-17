A = [1 -1.5 ; 0 1]; B = [0 ; 1];
SYS = ss(A,B,[],[]);

K1 = place(A,B,[-1+2i,-1-2i]);
K2 = place(A,B,[-1,-2]);
K3 = place(A,B,[1,2]);

sampling_intervals = 0.001:0.001:0.6';

max_poles = zeros(length(sampling_intervals),3);

for idx=1:length(sampling_intervals)
    h = sampling_intervals(idx);
    SYS_discrete = c2d(SYS,h);
    A_d = SYS_discrete.A;
    B_d = SYS_discrete.B;
    
    % maximum magnitude of discrete-time closed-loop poles for each controller
    max_pole_1 = max(abs(eig(A_d - B_d*K1))); 
    max_pole_2 = max(abs(eig(A_d - B_d*K2)));
    max_pole_3 = max(abs(eig(A_d - B_d*K3)));

    max_poles(idx,:) = [max_pole_1,max_pole_2,max_pole_3];
end

figure(1);
set(gca,'DefaultLineLineWidth',2)
set(gca,'Fontsize',12)
yline(1,'--');hold on;
plot(sampling_intervals,max_poles(:,1));hold on;
plot(sampling_intervals,max_poles(:,2));hold on;
plot(sampling_intervals,max_poles(:,3));hold on;
legend("","K1","K2","K3");
xlabel("sampling interval");
ylabel("magnitude");
title("Max. magnitude of closed-loop poles for various sampling intervals");
