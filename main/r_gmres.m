function out = r_gmres(A, b, param)

%% Function Definition %%
% R_GMRES - A function which solves a linear system Ax = b using the 
% (projected) recycled Generalized Minimum Residual Method (rGMRES).

%% Function Input %% 
%    A  : Coefficient matrix of linear system A x = b
%    b  : Right hand side of lienar system A x = b
%    param : A, struct with the following fields 
%          - n, number of rows and columns of A
%          - m, Arnoldi cycle length
%          - k, recycling subspace dimension
%          - U, matrix with columns forming basis for augmentation subspace
%          - C, matrix such that C = A*U;
%          - tol, Convergence tolerance of solver

%% Function Output %%
%    out : A struct with the following fields
%    - approx_sol, A vector containing a solution approximation for A x = b
%    - residuals, A vector containing the residual norm at the end of each
%                 Arnoldi cycle
%    - mv : Number of vectors matrix A is applied to in order to reach
%           convergence
%    - U (updated) matrix with columns forming basis for augmentation subspace
%    - C, matrix such that C = A*U;

n = param.n;  % Dimension of matrix A
m = param.m;  % Arnoldi cycle length
k = param.k;  % Dimension of Recycling Subspace

V = zeros(n,m+1);  % Stores Arnoldi basis
H = zeros(m+1,m);  % Stores Hessenberg matrix

% Stores relative residual norms at end of each cycle
residuals = zeros(1,param.max_cycles+1);  

% Right hand side for small linear system
small_rhs = zeros(m+1,1);  

% Vector to store solution approximation for A x = b
approx_sol = zeros(n,1);  

res = b;  % Initial residual vector (assuming initial null approximation)     
normr0 = norm(res); % Initial residual norm
normr = normr0;   
 
iter = 1;           % Iteration number
residnorm = 1;      % Initial relative residual norm is always 1

% Store initial relative residual in first entry of residuals vector
residuals(iter)= residnorm;   

% Initialize counter storing matrix vector products with A
mv = 0; 

% Ensure dimension of recycling subspace does not exceed Arnoldi 
% cycle length
if k > m
    warning('Recycling subspace must not exceed basis size.')
    k = m;
end

% If we do not yet have an initial recycling subspace, then for the first
% cycle we call the standard GMRES algorithm, which allows us to construct a 
% U for subsequent cycles
if isempty(param.U)
  
% mgs Arnoldi
V(:,1) = res/normr;

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

 % Solve small GMRES problem
 small_rhs(1) = normr;
 y = H\small_rhs;

 % Update solution and residual
 approx_sol =  V(:,1:m)*y;
 res = b - V*H*y;

 normr = norm(res);
 residnorm = normr/normr0;
 iter = iter+1;
 residuals(iter) = residnorm;
    
  % Construct recycling subspace for next cycle using harmonic Ritz vectors
  P = zeros(m,k);
  [harmVecs, harmVals] = eig(H'*H,H(1:m,1:m));
  harmVals = diag(harmVals);
    
  [~,iperm] = sort(abs(harmVals),'descend');
  for i = 1:k
      P(:,i) = harmVecs(:,iperm(i));
  end
    
  U = V(:,1:m) * P;
  C = A*U;

  % Impose columns of C be orthonormal
  [C,RR] = qr(C,0);
  U = U/RR;
 
  param.U = U;
  param.C = C;
    
  % Constructing C required an additional k applications of A
  mv = mv + k;

else


    U = param.U;
    C = param.C;

end

% Matrix to stores additional orthogonalization coefficients from 
% performing an Arnoldi cycle with a projected matrix
ortho_coeff = zeros(k, m);

while(residnorm > param.tol)
   
    tmp = C' * res;

    % Residual and solution projections
    r = res - C * tmp;
    approx_sol = approx_sol + U * tmp;

    normr = norm(r);
    V(:,1) = r/normr;

    % mgs Arnoldi 
     for j = 1:m
        v = A*V(:,j);
        mv = mv + 1;

        ortho_coeff(:,j) = C'*v;
        for i = 1:k
           v = v - C(:,i)*ortho_coeff(i,j);
        end
        
        for t = 1:j
            H(t,j) = V(:,t)'*v;
            v = v - V(:,t)*H(t,j);
        end

        [V(:,j+1),H(j+1,j)] = qr(v,0);
     end

    % Solve small GMRES problem
    small_rhs(1) = normr;
    y = H\small_rhs;


    % update solution and residual
    approx_sol = approx_sol + V(:,1:m) * y - U * ortho_coeff * y;
    res = r - V(:,1:m+1) * H * y;

    normr = norm(res);
    residnorm = normr/normr0;

    iter = iter + 1;
    residuals(iter) = residnorm;

    % Update recycle space using harmonic ritz vectors
    [U,C] = r_gmres_harmritz_recycling(A,param,ortho_coeff,V,H);
    param.U = U;
    param.C = C;
    mv = mv + k;

end

    % output quantities into out struct
    out.U = param.U;
    out.C = param.C;
    out.approx_sol = approx_sol;
    out.residuals = residuals;
    out.mv = mv;

end
