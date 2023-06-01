K = [-16/3 4 0];
h_values = 0.01 : 0.01 : 0.6;
stability_bools_tohold = arrayfun(@(h) is_system_stable_tohold(h,K), h_values);
stability_bools_tozero = arrayfun(@(h) is_system_stable_tozero(h,K), h_values);

figure();
subplot(2,1,1); scatter(h_values,stability_bools_tohold); ylim([-0.1,1.1]); title("To hold"); xlabel("sampling interval (h)"); grid on;
subplot(2,1,2); scatter(h_values,stability_bools_tozero); ylim([-0.1,1.1]); title("To zero"); xlabel("sampling interval (h)"); grid on;

function stab = is_system_stable_tohold(h,K)
    h
    A = [1 -1.5 ; 0 1]; B = [0 ; 1];
    A_d = expm(A*h); B_d = (expm(A*h) - eye(2))*(A\B);
    F_0 = [A_d zeros(2,1) ; zeros(1,2) 0]; G_0 = [B_d ; 1];
    F_1 = [A_d B_d ;  zeros(1,2) 1];       G_1 = [zeros(2,1) ; 0];
    M_0 = F_0 - G_0*K; M_1 = F_1 - G_1*K;

    M_000 = M_0^3;
    M_001 = M_1 * M_0^2;

    P = sdpvar(3,3);
    constraints = [P >= 0, ...
                   M_000' * P * M_000 - P <= 0, ...
                   M_001' * P * M_001   - P <= 0];
    options = sdpsettings('verbose',0);
    diagnostics = optimize(constraints,P,options);
    stab = strcmp('Successfully solved (LMILAB)',diagnostics.info);
end

function stab = is_system_stable_tozero(h,K)
    h
    A = [1 -1.5 ; 0 1]; B = [0 ; 1];
    A_d = expm(A*h); B_d = (expm(A*h) - eye(2))*(A\B);
    F_0 = [A_d zeros(2,1) ; zeros(1,2) 0]; G_0 = [B_d ; 1];
    F_1 = [A_d zeros(2,1) ; zeros(1,2) 1]; G_1 = [zeros(2,1) ; 0];
    M_0 = F_0 - G_0*K; M_1 = F_1 - G_1*K;

    M_000 = M_0^3;
    M_001 = M_1 * M_0^2;

    P = sdpvar(3,3);
    constraints = [P >= eye(3), ...
                   M_000' * P * M_000 - P <= 0, ...
                   M_001' * P * M_001   - P <= 0];
    options = sdpsettings('verbose',0);
    diagnostics = optimize(constraints,P,options);
    stab = strcmp('Successfully solved (LMILAB)',diagnostics.info);
end