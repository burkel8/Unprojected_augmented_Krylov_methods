%%%%%         gmres_residual_approx_test.m    %%%%%
%   Computes the exact residual computation in unprojected recycled GMRES (ur_gmres) 
%   and compares it to the residual estimate proposed in the paper.

%   The test matrix is a QCD matrix of size 3072 x 3072 

%%%%% User defined parameters to be tuned are defined here  %%%

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

rGMRESp = p;
ur_gmres_p = p;

fprintf("\n Solving a linear system using urGMRES \n");
pause(5);

% Create new random rhs for each system    
rng(4);
b = randn(n,1);

% Call ur_gmres on each system
ur_gmres_o = ur_gmres(A, b, ur_gmres_p);
ur_gmres_p.U = ur_gmres_o.U;
ur_gmres_p.C = ur_gmres_o.C;

fprintf("\n Plotting estimated residual norm vs exact residual norm for last system \n");

% Plot residual estimate vs true residual
semilogy(ur_gmres_o.residuals_approx,'LineWidth',3);
hold on; 
semilogy(ur_gmres_o.residuals,'LineWidth',3);
hold off;
legend('Estimated residual norm','Exact residual norm','FontSize',20);
xlabel("Restart Number");
set(gca,"FontSize",15)
grid on;
