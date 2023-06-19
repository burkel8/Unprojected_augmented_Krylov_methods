%%%%%         gmres_residual_approx_test.m    %%%%%
%   Computes the exact residual computation in unprojected recycled GMRES (ur_gmres) 
%   and compares it to the residual estimate proposed in the paper.

%   The test matrix is a QCD matrix of size 3072 x 3072 

%%%%% User defined parameters to be tuned are defined here  %%%

% p is a struct with various fields

p.m = 30;           % Dimension of Krylov subspace
p.max_cycles = 5;   % Max number of Arnoldi cycles
p.k = 10;           % Recycling subspace dimension
p.tol = 1e-15;      % Convergence Tolerance
num_systems = 4;    % Number of linear systems in a sequence
p.U = [];       % Recycling subspace basis
p.C = [];       % C such that C = A*U;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QCD matrix
addpath(genpath('../'));
load("smallLQCD_A1.mat");
n = size(A1,1);
A = A1 - 0.65*speye(n);
p.n = n;

rGMRESp = p;
ur_gmres_p = p;

% Solve all systems
for i = 1:num_systems

% Create new random rhs for each system    
rng(i);
b = randn(n,1);

% Call ur_gmres on each system
ur_gmres_o = ur_gmres(A, b, ur_gmres_p);
ur_gmres_p.U = ur_gmres_o.U;
ur_gmres_p.C = ur_gmres_o.C;

end

% Plot residual estimate vs true residual
semilogy(ur_gmres_o.residuals_approx,'LineWidth',2);
hold on; 
semilogy(ur_gmres_o.residuals,'LineWidth',2);
hold off;

legend('Unprojected rGMRES residual estimate','Unprojected rGMRES residual','FontSize',12);
xlabel("Restart Number");
ylabel("Residuals");
grid on;
