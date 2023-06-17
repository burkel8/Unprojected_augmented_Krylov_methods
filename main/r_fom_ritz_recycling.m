function [U,C] = r_fom_ritz_recycling(A,param,ortho_coeff,V,H)

%% Function Definition %%
% R_FOM_RITZ_RECYCLING - A function which construct a recycling
% subspace for rFOM using an augmented Krylov subspace K + U, by solving
% the appropriate Ritz problem for rFOM.

%% Function Input %% 
%    A  : Coefficient matrix of linear system A x = b
%    param : A, struct with the following fields 
%          - n, number of rows and columns of A
%          - m, Arnoldi cycle length
%          - k, recycling subspace dimension
%          - U, matrix with columns forming basis for augmentation subspace
%          - C, matrix such that C = A*U;
%          - tol, Convergence tolerance of solver
%       L  : matrix storing orthogonalization coefficients from additional
%            orthogonalizations in Arnoldi
%       V   : matrix whos columns form a basis for Krylov subspace K
%       H   : Hessenberg matrix from Arnoldi to build basis for K
%       

%% Function Output %%
%    - U (updated) matrix with columns forming basis for augmentation subspace
%    - C, matrix such that C = A*U;

    U = param.U;
    C = param.C;
    m = param.m;
    k = param.k;
    P = zeros(m + k, k);
    
    Vhat = [U V(:,1:m)];
    UTC = U'*C;

    % Set up and solve harmonic Ritz problem
    A1 = eye(k+m);

    A2 = zeros(k + m,k + m);
    A2(1:k,1:k) = UTC;
    A2(1:k,k+1:end) = UTC*ortho_coeff;
    A2(k+1:end,k+1:end) = H(1:m,1:m);
    
    [harmVecs, harmVals] = eig(A1,A2);
    harmVals = diag(harmVals);
    
    [~,iperm] = sort(abs(harmVals),'descend');
    for i = 1:k
        P(:,i) = harmVecs(:,iperm(i));
    end
    
    U = Vhat * P;

    % Impose orthonormal columns of U
    [U,~]=qr(U,0);
    C = A*U;

end