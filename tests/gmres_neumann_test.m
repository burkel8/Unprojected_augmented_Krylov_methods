%%%%%         gmres_neumann_test.m    %%%%%
%   Tests and compares unprojected recycled GMRES (ur_gmres) to 
%   recycled GMRES (r_gmres) and standard GMRES (gmres). Convergence curves
%   and number of vectors the coefficient matrix A is applied to is recorded 
%   to compare each algorithm.

%   The test matrix is a Neumann matrix of size 22500 x 22500 

%%%%% User defined parameters to be tuned are defined here  %%%

% p is a struct with various fields
p.m = 90;           % Dimension of Krylov subspace
p.max_cycles = 5;   % Max number of Arnoldi cycles
p.k = 20;           % Reycling subspace dimension
p.tol = 1e-14;      % Convergence Tolerance
num_systems = 5;    % Number of linear systems in a sequence
p.U = [];       % Recycling subspace basis
p.C = [];       % C such that C = A*U;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct Neumann Matrix
addpath(genpath('../'));
rng(4);
n = 22500;
A = gallery('neumann', n) + 0.0001*speye(n);

% vectors to store number of vectors matrix A is applied to for each
% problem in the sequence
gmres_mv = zeros(1,num_systems);
r_gmres_mv = zeros(1,num_systems);
ur_gmres_mv = zeros(1,num_systems);

p.n = n;  %  A is n x n matrix

% Input struct for each method will be the same
gmres_p = p;
r_gmres_p = p;
ur_gmres_p = p;

% Counters to record total number of vectors A is applied to for each
% method
tot_gmres_mv = 0;
tot_r_gmres_mv = 0;
tot_ur_gmres_mv = 0;

fprintf("\n Solving a sequence of %d linear system(s) using GMRES, rGMRES and urGMRES\n", num_systems);
fprintf("\n  Printing Number of MATVEC's required for each method to converge \n")
pause(5);

% Loop through the full sequence of systems and solve each using the three 
% methods (gmres , r_gmres and ur_gmres)
for i = 1:num_systems

fprintf("\n #######  System %d #######  \n", i);

% Create random right hand size for each system
rng(i);
b = randn(n,1);

%% Call gmres
fprintf("\n Running GMRES  \n");
gmres_o = gmres(A, b, gmres_p);

%% Call r_gmres
fprintf("\n Running rGMRES  \n");
r_gmres_o = r_gmres(A, b, r_gmres_p);

% output recycling subspace from r_gmres call and 
% add it as an input for next system call to r_gmres

r_gmres_p.U = r_gmres_o.U;
r_gmres_p.C = r_gmres_o.C;

%% Call ur_gmres
fprintf("\n Running urGMRES  \n");
ur_gmres_o = ur_gmres(A, b, ur_gmres_p);

% output recycling subspace from ur_gmres call and 
% add it as an input for next system call to ur_gmres

ur_gmres_p.U = ur_gmres_o.U;
ur_gmres_p.C = ur_gmres_o.C;

% Accumalate number of A applications for all methods
tot_gmres_mv = tot_gmres_mv + gmres_o.mv;
tot_r_gmres_mv = tot_r_gmres_mv + r_gmres_o.mv;
tot_ur_gmres_mv = tot_ur_gmres_mv + ur_gmres_o.mv;

gmres_mv(1,i) = tot_gmres_mv;
r_gmres_mv(1,i) = tot_r_gmres_mv;
ur_gmres_mv(1,i) = tot_ur_gmres_mv;

fprintf("\n             MATVEC's            \n");
fprintf('\n GMRES %d rGMRES %d urGMRES %d \n',gmres_o.mv,r_gmres_o.mv, ur_gmres_o.mv);
pause(5);

end

fprintf("\n ######## Total MATVEC's #######  \n");
fprintf("\n GMRES %d rGMRES %d urGMRES %d\n", tot_gmres_mv,tot_r_gmres_mv, tot_ur_gmres_mv);

% plot convergence curve of final system.
semilogy(gmres_o.residuals,'LineWidth',4);
hold on;
semilogy(r_gmres_o.residuals,'LineWidth',4);
hold on; 
semilogy(ur_gmres_o.residuals,'LineWidth',4);
hold off;
legend('GMRES','rGMRES','urGMRES','FontSize',20);
xlabel("Restart Number");
ylabel("Relative Residual");
set(gca,"FontSize",15)
grid on;

