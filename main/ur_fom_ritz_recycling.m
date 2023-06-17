function [U,C] = ur_fom_ritz_recycling(A,param,V,H)

%% Function Definition %%
% UR_FOM_HARMRITZ_RECYCLING - A function which construct a recycling
% subspace for urFOM using an augmented Krylov subspace K + U, by solving
% the appropriate Ritz problem for urFOM.

%% Function Input %% 
%    A  : Coefficient matrix of linear system A x = b
%    param : A, struct with the following fields 
%          - n, number of rows and columns of A
%          - m, Arnoldi cycle length
%          - k, recycling subspace dimension
%          - U, matrix with columns forming basis for augmentation subspace
%          - C, matrix such that C = A*U;
%          - tol, Convergence tolerance of solver
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
    UTW = U'*V(:,1:m);
    WTU = V(:,1:m)'*U;
    WTC = V(:,1:m)'*C;

    % Set up and solve Ritz problem
    A1 = eye(k+m);
    A1(1:k,1:k) = U'*U; 
    A1(1:k,k+1:end) = UTW;
    A1(k+1:end,1:k) = WTU;

    A2 = zeros(k + m,k + m);
    A2(1:k,1:k) = UTC;
    A2(1:k,k+1:end) = U'*V*H;
    A2(k+1:end,1:k) = WTC;
    A2(k+1:end,k+1:end) = H(1:m,1:m);
    
    [harmVecs, harmVals] = eig(A1,A2);
    harmVals = diag(harmVals);
    
    [~,iperm] = sort(abs(harmVals),'descend');
    for i = 1:k
        P(:,i) = harmVecs(:,iperm(i));
    end
    
    U = Vhat * P;
    C = A*U;

    % Impose C to have orthonormal columns to monitor residual norm cheaply
    [C,RR] = qr(C,0);
    U = U/RR;

end