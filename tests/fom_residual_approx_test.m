%%%%%         fom_residual_approx_test.m    %%%%%
%   Computes the exact residual computation in unprojected recycled GMRES (ur_fom) 
%   And compares it to the residual estimate proposed in the paper.

%   The test matrix is a QCD matrix of size 3072 x 3072 

% p is a struct with various fields
p.m = 120;           % Dimension of Krylov subspace
p.max_cycles = 5;   % Max number of Arnoldi cycles
p.k = 20;           % Recycling subspace dimension
p.tol = 1e-15;      % Convergence Tolerance
num_systems = 1;    % Number of linear systems in a sequence
p.U = [];       % Recycling subspace basis
p.C = [];       % C such that C = A*U;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 120;
n = N*N;
A = gallery("poisson",N);
p.n = n;

r_fom_mv = zeros(1,num_systems);
ur_fom_mv = zeros(1,num_systems);

p.n = n;

r_domp = p;
ur_fom_p = p;

fprintf("\n Solving a linear system using urFOM\n");
pause(5);

% Create new rhs for each system
rng(4);
b = randn(n,1);

% Call ur_fom for each system
ur_fom_o = ur_fom(A, b, ur_fom_p);
ur_fom_p.U = ur_fom_o.U;
ur_fom_p.C = ur_fom_o.C;

% Plot residual estimate vs true residual
fprintf("\n Plotting estimated residual norm vs exact residual norm for last system \n");
semilogy(ur_fom_o.residuals_approx,'LineWidth',2);
hold on; 
semilogy(ur_fom_o.residuals,'LineWidth',2);
hold off;

legend('urFOM estimated residual norm','urFOM exact residual norm','FontSize',12);
xlabel("Restart Number");
ylabel("Residuals");
grid on;
