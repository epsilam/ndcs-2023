% Goal: find a controller K such that the possible selection of h_n
% which stabilizes both systems with 
% (h,tau)=(0.5h_n,0.25h_n) and (h,tau)=(1.5h_n,0.75h_n)
% is maximized

%[K,fval,exitflag,output,population,scores] = ga(@objective_function,3);
%population
%K
%fval

%[K,fval] = fminsearch(@objective_function,[-1.3531    3.9921    0.9982]);

test_stability(0.5344,[-1.3531    3.9921    0.9982])

function out = objective_function(K)
    h_n_previous = 0.01;
    for h_n=0.01:0.01:2
        stab = test_stability(h_n,K);
        if ~stab
            out = -h_n_previous;
            break
        end
        h_n_previous = h_n;
    end
end

function stability = test_stability(h_n,K)
    A = [1 -1.5 ; 0 1]; B = [0 ; 1];

    h   = 0.5*h_n;
    tau = 0.25*h_n;
    F_x = expm(A*h); F_u = (expm(A*h) - expm(A*(h-tau)))*(A\B);
    G_1 = (expm(A*(h-tau)) - eye(2))*(A\B);
    A_aug = [F_x F_u ; zeros(1,2) 0 ];
    B_aug = [G_1 ; 1];
    stability_1 = max(abs(eig(A_aug - B_aug*K))) < 1;

    % we require both systems to be 
    if stability_1 
        h   = 1.5*h_n;
        tau = 0.75*h_n;
        F_x = expm(A*h); F_u = (expm(A*h) - expm(A*(h-tau)))*(A\B);
        G_1 = (expm(A*(h-tau)) - eye(2))*(A\B);
        A_aug = [F_x F_u ; zeros(1,2) 0 ];
        B_aug = [G_1 ; 1];
        stability_2 = max(abs(eig(A_aug - B_aug*K))) < 1;
        stability = stability_2;
    else 
        stability = 0;
    end
end