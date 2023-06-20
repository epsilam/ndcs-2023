aircraft;



zero = zeros(n,m*Tfinal);

% Matrix and vector for linear equality constraint


% Find a point in the constraint set, which is necessarily in its
% relative interior
u = -pinv(A_h)*b_h;

u_test_lower = -12 < u;
u_test_upper = u < 12;

if ~nnz(u_test_lower==0) && ~nnz(u_test_upper==0)
    disp("Slater's condition holds")
end

function [A_h, b_h] equality_constraints(A1,A2,A3,A4,B1,B2,B3,B4,x01,x02,x03,x04,Tfinal)
    C1 = Ci(A1,B1,Tfinal); v1 = vi(A1,x01,Tfinal);
    C2 = Ci(A2,B2,Tfinal); v2 = vi(A2,x02,Tfinal);
    C3 = Ci(A3,B3,Tfinal); v3 = vi(A3,x03,Tfinal);
    C4 = Ci(A4,B4,Tfinal); v4 = vi(A4,x04,Tfinal);
    
    A_h = [C1,   -C2,  zero, zero ;
     zero,  C2,  -C3,  zero ;
     zero, zero, C3,   C4   ];
    b_h = [v1 - v2 ; v2 - v3 ; v3 - v4];
end

function out = Ci(Ai,Bi,Tfinal)
    n = size(Ai,1); m = size(Bi,2);
    out = zeros(n,m*Tfinal);
    for t=1:Tfinal
        out(:, m*t+1:m*(t+1)) = Ai^(Tfinal-t)*Bi;
    end
end

function out = vi(Ai,x0i,Tfinal)
    out = Ai^Tfinal * x0i;
end