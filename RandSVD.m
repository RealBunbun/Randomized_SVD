function [U, S, V] = RandSVD(A, Rank)

N_Oversamples    = 10; 
N_Subspace_Iters = 5;
Return_Range     = 'False';

if N_Oversamples ==0
  N_Samples = 2*Rank;
else
  N_Samples = Rank + N_Oversamples;
end

% // Phase 1
Q = Find_Range(A, N_Samples, N_Subspace_Iters);

% // Phase 2
B = transpose(Q) * A;
[U_tilde, S, V] = svd(B, 'econ');
U = Q * U_tilde;

% // Truncation
U = U(:,1:Rank);
S = diag(S);
S = S(1:Rank);
V = V(:,1:Rank);

end

function [Q] = Find_Range(A, N_Samples, N_Subspace_Iters)
[~,n] = size(A);
O = randn(n, N_Samples);
Y = A * O;

if N_Subspace_Iters==0
  Q = Ortho_Basis(Y);
else
  Q = Subspace_Iter(A, Y, N_Subspace_Iters);
end

end

function [Q] = Ortho_Basis(Y)
[Q,~] = qr(Y,0);
end

function [Q] = Subspace_Iter(A, Y0, N_Iters)

Q = Ortho_Basis(Y0);
for ii = 1:N_Iters
  Z = Ortho_Basis(transpose(A) * Q);
  Q = Ortho_Basis(A * Z);
end

end

