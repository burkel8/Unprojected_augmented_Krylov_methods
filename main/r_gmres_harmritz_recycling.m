function [U,C] = r_gmres_harmritz_recycling(A,param,L,V,H)

%% Function Definition %%
% R_GMRES_HARMRITZ_RECYCLING - A function which construct a recycling
% subspace for rGMRES using an augmented Krylov subspace K + U, by solving
% the appropriate harmonic Ritz problem for rGMRES.

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
    CTU = C'*U;

    % Set up and solve appropriate harmonic Ritz problem
    A1 = zeros(k+m);
    A1(1:k,1:k) = CTU;
    A1(k+1:end,1:k) = (V*H + C*L)'*U;
    A1(k+1:end,k+1:end) = H(1:m,1:m)';

    A2 = eye(k+m);
    A2(1:k,k+1:end) = L;
    A2(k+1:end,1:k) = L';
    A2(k+1:end,k+1:end) = H'*H + L'*L;
    
    [harmVecs, harmVals] = eig(A1,A2);
    harmVals = diag(harmVals);
    
    % Select k "smallest" eigenvectors
    [~,iperm] = sort(abs(harmVals),'descend');
    for i = 1:k
        P(:,i) = harmVecs(:,iperm(i));
    end
    
    U = Vhat * P;
    C = A*U;

    % Impose orthonormal columns of C for rGMRES
    [C,RR] = qr(C,0);
    U = U/RR;

end