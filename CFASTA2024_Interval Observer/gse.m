function [K,V,W,Nhat,Dhat] = gse(A,B,F,Z)
%GSE solves Generalized Sylvester Equations (GSEs) AV + BW = VF and 
%produces state feedback gain K for linear time-invariant systems.
%   **INPUT**
%               A: system matrix
%               B: input matrix
%               s: prescribed poles
%               F: matrix with expected poles on the diagonal
%               Z: free matrix
%   **OUTPUT**
%               K: state feedback gain matrix
%               V: solution to GSEs
%               W: solution to GSEs
%               Nhat: solution to (sI-A)Nhat = B Dhat
%               Dhat: solution to (sI-A)Nhat = B Dhat
%   Consider linear dynamic systems with state x,
%           continuous case: dx/dt = Ax + Bu
%           discrete case:   x[k+1] = Ax[k] + Bu[k]
%   The control input is designed as u = Kx (u[k] = Kx[k]), where K = W / V.

% -------------------------------------------------------------------------
% Version:              1.0
% Author:               Weijie Ren
% Contact:              weijie.ren@outlook.com
% Initial modified:     Nov. 14, 2023
% Last modified:        
% -------------------------------------------------------------------------

% Dimension
n = size(A,1);
r = size(B,2);
% Declare array
P = cell(n,1);
Q = cell(n,1);
Nhat = cell(n,1);
Dhat = cell(n,1);
V = zeros(n,n);
W = zeros(r,n);
% Compute Nhat and Dhat
s = diag(F);
for i = 1:n
    % Singular Value Decomposition [sI-A  -B]
    [UU,~,VV] = svd([s(i)*eye(n)-A, -B]);
    P{i} = UU';
    Q{i} = VV;
    % Split Q
    Nhat{i} = Q{i}(1:n,n+1:n+r);
    Dhat{i} = Q{i}(n+1:n+r,n+1:n+r);
    % Compute V and W
    V(:,i) = Nhat{i}*Z(:,i);
    W(:,i) = Dhat{i}*Z(:,i);
end
% Calculate control gain
K = W / V;
end