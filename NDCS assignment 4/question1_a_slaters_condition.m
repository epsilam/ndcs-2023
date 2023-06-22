% Import standard function definitions
question1_a_function_definitions;

% Find a point in the constraint set, which is necessarily in its
% relative interior
u = -pinv(A_h)*b_h;

u_test_lower = -12 < u;
u_test_upper = u < 12;

if ~nnz(u_test_lower==0) && ~nnz(u_test_upper==0)
    disp("Slater's condition holds.")
else
    disp("Slater's condition does not hold.")
end