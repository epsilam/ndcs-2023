K = [-16/3 4 0];
h_values = 0.05 : 0.01 : 0.5;
stability_bools_tohold = arrayfun(@(h) is_system_stable_tohold(h,K), h_values);
stability_bools_tozero = arrayfun(@(h) is_system_stable_tozero(h,K), h_values);

figure();
subplot(2,1,1); scatter(h_values,stability_bools_tohold,"filled"); ylim([-0.1,1.1]); title("To hold"); xlabel("sampling interval (h)"); grid on;
subplot(2,1,2); scatter(h_values,stability_bools_tozero,"filled"); ylim([-0.1,1.1]); title("To zero"); xlabel("sampling interval (h)"); grid on;

function stab = is_system_stable_tohold(h,K)
    h
    A = [1 -1.5 ; 0 1]; B = [0 ; 1];
    A_d = expm(A*h); B_d = (expm(A*h) - eye(2))*(A\B);
    F_0 = [A_d zeros(2,1) ; zeros(1,2) 0]; G_0 = [B_d ; 1];
    F_1 = [A_d B_d ;  zeros(1,2) 1];       G_1 = [zeros(2,1) ; 0];
    M_0 = F_0 - G_0*K; M_1 = F_1 - G_1*K;
    stab = solve_LMIs(M_0,M_1);
end

function stab = is_system_stable_tozero(h,K)
    h
    A = [1 -1.5 ; 0 1]; B = [0 ; 1];
    A_d = expm(A*h); B_d = (expm(A*h) - eye(2))*(A\B);
    F_0 = [A_d zeros(2,1) ; zeros(1,2) 0]; G_0 = [B_d ; 1];
    F_1 = [A_d zeros(2,1) ; zeros(1,2) 1]; G_1 = [zeros(2,1) ; 0];
    M_0 = F_0 - G_0*K; M_1 = F_1 - G_1*K;
    stab = solve_LMIs(M_0,M_1);
end

function stab = solve_LMIs(M_0,M_1)
    P_0xx = sdpvar(3,3);
    P_1xx = sdpvar(3,3);
    P_x0x = sdpvar(3,3);
    P_x1x = sdpvar(3,3);
    P_xx0 = sdpvar(3,3);
    P_xx1 = sdpvar(3,3);

    constraints = [P_0xx >= eye(3), ...
                   P_1xx >= eye(3), ...
                   P_x0x >= eye(3), ...
                   P_x1x >= eye(3), ...
                   P_xx0 >= eye(3), ...
                   P_xx1 >= eye(3), ...
                   P_0xx - M_0' * (0.99*P_xx0 + 0.99*P_xx1) * M_0 >= eye(3), ...
                   P_1xx - M_1' * (0.01*P_xx0 + 0.01*P_xx1) * M_1 >= eye(3), ...
                   P_x0x - M_0' * (0.99*P_0xx + 0.99*P_1xx) * M_0 >= eye(3), ...
                   P_x1x - M_1' * (0.01*P_0xx + 0.01*P_1xx) * M_1 >= eye(3), ...
                   P_xx0 - M_0' * (0.49*P_x0x + 0.49*P_x1x) * M_0 >= eye(3), ...
                   P_xx1 - M_1' * (0.51*P_x0x + 0.51*P_x1x) * M_1 >= eye(3)];
    options = sdpsettings('verbose',0);
    diagnostics = optimize(constraints,[P_0xx,P_1xx,P_x0x,P_x1x,P_xx0,P_xx1],options);
    %diagnostics = optimize(constraints,[]);
    %disp("Feasibility:::")
    %disp(diagnostics.info)
    stab = strcmp('Successfully solved (LMILAB)',diagnostics.info);
end