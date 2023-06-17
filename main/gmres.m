function out = gmres(A, b, param)

%% Function Definition %%
% GMRES - A function which solves a linear system Ax = b using the standard 
% restarted Generalized Minimum Residual Method (GMRES).

%% Function Input %% 
%    A  : Coefficient matrix of linear system A x = b
%    b  : Right hand side of lienar system A x = b
%    param : A, struct with the following fields 
%          - n, number of rows and columns of A
%          - m, Arnoldi cycle length
%          - tol, Convergence tolerance of solver

%% Function Output %%
%    out : A struct with the following fields
%    - approx_sol, A vector containing a solution approximation for A x = b
%    - residuals, A vector containing the residual norm at the end of each
%                 Arnoldi cycle
%    - mv : Number of vectors matrix A is applied to in order to reach
%           convergence

n = param.n;  % Dimension of matrix A
m = param.m;  % Arnoldi cycle length

V = zeros(n,m+1);  % Stores Arnoldi basis
H = zeros(m+1,m);  % Stores Hessenberg matrix

% Stores relative residual norms at end of each cycle
residuals = zeros(1,param.max_cycles+1);  

% Right hand side for small linear system
small_rhs = zeros(m+1,1);  

% Vector to store solution approximation for A x = b
approx_sol = zeros(n,1);  

r = b;  % Initial residual vector (assuming initial null approximation)     
normr0 = norm(r); % Initial residual norm
normr = normr0;   
 
iter = 1;           % Iteration number
residnorm = 1;      % Initial relative residual norm is always 1

% Store initial relative residual in first entry of residuals vector
residuals(iter)= residnorm;   

% Initialize counter storing matrix vector products with A
mv = 0; 

% Run restarted FOM until convergence

while(residnorm > param.tol)

V(:,1) = r/normr;

% mgs Arnoldi
for j=1:m

   V(:,j+1) = A*V(:,j);
   mv = mv + 1;
   
   for i=1:j
       H(i,j)= V(:,i)'*V(:,j+1);
       V(:,j+1)= V(:,j+1) - V(:,i)*H(i,j);
   end
    
    H(j+1,j) = norm(V(:,j+1));
    V(:,j+1) = V(:,j+1)/H(j+1,j);
 end

    % Solve small FOM linear system
    small_rhs(1) = normr;
    y = H\small_rhs;

    % Update solution and residual
    approx_sol = approx_sol + V(:,1:m) * y;
    r = r - V(:,1:m+1) * H * y;
   
    normr = norm(r);
    residnorm = normr/normr0;

    iter = iter + 1;
    residuals(iter) = residnorm;
   
end
    
    % Output approx_sol and residuals vector in param struct
    out.approx_sol = approx_sol;
    out.residuals = residuals;
    out.mv = mv;
   
end
