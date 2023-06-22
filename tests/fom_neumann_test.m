%%%%%         fom_neumann_test.m    %%%%%
%   Tests and compares unprojected recycled GMRES (ur_fom) to 
%   recycled GMRES (r_fom) and standard FOM (fom). Convergence curves
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
fom_mv = zeros(1,num_systems);
r_fom_mv = zeros(1,num_systems);
ur_fom_mv = zeros(1,num_systems);

p.n = n;  %  A is n x n matrix

% Input struct for each method will be the same
fom_p = p;
r_fom_p = p;
ur_fom_p = p;

% Counters to record total number of vectors A is applied to for each
% method
tot_fom_mv = 0;
tot_r_fom_mv = 0;
tot_ur_fom_mv = 0;

fprintf("\n Solving a sequence of %d linear system(s) using FOM, rFOM and urFOM\n", num_systems);
fprintf("\n  Printing Number of MATVEC's required for each method to converge \n")
pause(5);

% Loop through the full sequence of systems and solve each using the three 
% methods (fom , r_fom and ur_fom)
for i = 1:num_systems

fprintf("\n #######  System %d #######  \n", i);

% Create random right hand size for each system
rng(i);
b = randn(n,1);

%% Call fom
fprintf("\n Running FOM  \n");
fom_o = fom(A, b, fom_p);

%% Call r_fom
fprintf("\n Running rFOM  \n");
r_fom_o = r_fom(A, b, r_fom_p);

% output recycling subspace from r_fom call and 
% add it as an input for next system call to r_fom

r_fom_p.U = r_fom_o.U;
r_fom_p.C = r_fom_o.C;

%% Call ur_fom
fprintf("\n Running urFOM  \n");
ur_fom_o = ur_fom(A, b, ur_fom_p);

% output recycling subspace from ur_fom call and 
% add it as an input for next system call to ur_fom

ur_fom_p.U = ur_fom_o.U;
ur_fom_p.C = ur_fom_o.C;

% Accumalate number of A applications for all methods
tot_fom_mv = tot_fom_mv + fom_o.mv;
tot_r_fom_mv = tot_r_fom_mv + r_fom_o.mv;
tot_ur_fom_mv = tot_ur_fom_mv + ur_fom_o.mv;

fom_mv(1,i) = tot_fom_mv;
r_fom_mv(1,i) = tot_r_fom_mv;
ur_fom_mv(1,i) = tot_ur_fom_mv;

fprintf("\n             MATVEC's            \n");
fprintf('\n FOM %d rFOM %d urFOM %d \n',fom_o.mv,r_fom_o.mv, ur_fom_o.mv);
pause(5);

end

fprintf("\n ######## Total MATVEC's #######  \n");
fprintf("\n FOM %d rFOM %d urFOM %d\n", tot_fom_mv,tot_r_fom_mv, tot_ur_fom_mv);

% plot convergence curve of final system.
semilogy(fom_o.residuals,'LineWidth',4);
hold on;
semilogy(r_fom_o.residuals,'LineWidth',4);
hold on; 
semilogy(ur_fom_o.residuals,'LineWidth',4);
hold off;
legend('FOM','rFOM','urFOM','FontSize',20);
xlabel("Restart Number");
ylabel("Relative Residual");
set(gca,"FontSize",15)
grid on;