h_avg = 0.3006; % largest inter-sample time obtained in question 5.2

h = h_avg;
A = [1 -1.5 ; 0 1]; B = [0 ; 1]; K = [-16/3 4];
A_d = expm(A*h); B_d = (expm(A*h) - eye(2))*(A\B);

max(abs(eig(A_d - B_d * K))) < 1

