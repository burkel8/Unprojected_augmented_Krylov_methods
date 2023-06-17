%%%%%         fom_residual_approx_test.m    %%%%%
%   Computes the exact residual computation in unprojected recycled GMRES (ur_fom) 
%   And compares it to the residual estimate proposed in the paper.

%   The test matrix is a QCD matrix of size 3072 x 3072 

%%%%% User defined parameters to be tuned are defined here  %%%
p.m = 30;           % Dimension of Krylov subspace
p.max_cycles = 5;   % Max number of Arnoldi cycles
p.k = 10;           % Recycling subspace dimension
p.tol = 1e-15;      % Convergence Tolerance
num_systems = 4;    % Number of linear systems in a sequence
p.U = [];       % Recycling subspace basis
p.C = [];       % C such that C = A*U;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QCD test
addpath(genpath('../'));
load("smallLQCD_A1.mat");
n = size(A1,1);
A = A1 - 0.65*speye(n);

r_fom_mv = zeros(1,num_systems);
ur_fom_mv = zeros(1,num_systems);

p.n = n;

r_domp = p;
ur_fom_p = p;

% Solve all systems
for i = 1:num_systems

% Create new rhs for each system
rng(i);
b = randn(n,1);

% Call ur_fom for each system
ur_fom_o = ur_fom(A, b, ur_fom_p);
ur_fom_p.U = ur_fom_o.U;
ur_fom_p.C = ur_fom_o.C;

end

% Plot residual estimate vs true residual
semilogy(ur_fom_o.residuals_approx,'LineWidth',2);
hold on; 
semilogy(ur_fom_o.residuals,'LineWidth',2);
hold off;

legend('Unprojected r_dom residual estimate','Unprojected r_dom residual','FontSize',12);
xlabel("Restart Number");
ylabel("Residuals");
grid on;
