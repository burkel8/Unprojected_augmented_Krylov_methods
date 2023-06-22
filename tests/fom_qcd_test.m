%%%%%         fom_qcd_test.m    %%%%%
%   Tests and compares unprojected recycled FOM (ur_fom) to 
%   recycled FOM (r_fom) and standard FOM (fom). Convergence curves
%   and number of vectors the coefficient matrix A is applied to is recorded 
%   to compare each algorithm.

%   The test matrix is a QCD matrix of size 3072 x 3072 

%%%%% User defined parameters to be tuned are defined here  %%%

% p is a struct with various fields
p.m = 50;           % Dimension of Krylov subspace
p.max_cycles = 5;   % Max number of Arnoldi cycles
p.k = 10;           % Recycling subspace dimension
p.tol = 1e-13;      % Convergence Tolerance
num_systems = 3;    % Number of linear systems in a sequence
p.U = [];       % Recycling subspace basis
p.C = [];       % C such that C = A*U;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../'))
load("smallLQCD_A1.mat");
n = size(A1,1);
A = A1 - 0.65*speye(n);

% Counters to record total number of vectors A is applied to for each
% method
fom_mv = zeros(1,num_systems);
r_fom_mv = zeros(1,num_systems);
ur_fom_mv = zeros(1,num_systems);

p.n = n; %  A is n x n matrix

% Input struct for each method will be the same
fom_p = p;
r_fom_p = p;
ur_fom_p = p;

%Record tot number of A applications for both methods
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

fprintf("\n Running FOM  \n");
fom_o = fom(A, b, fom_p);

fprintf("\n Running rFOM  \n");
r_fom_o = r_fom(A, b, r_fom_p);

% output recycling subspace from r_fom call and 
% add it as an input for next system call to r_fom

r_fom_p.U = r_fom_o.U;
r_fom_p.C = r_fom_o.C;

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

A = A + 0.0001*sprand(A);

end

fprintf("\n ######## Total MATVEC's #######  \n");
fprintf("\n FOM %d rFOM %d urFOM %d\n", tot_fom_mv,tot_r_fom_mv, tot_ur_fom_mv);

% plot convergence curve of final system.

h = semilogy(fom_o.residuals,'LineWidth',4);
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
