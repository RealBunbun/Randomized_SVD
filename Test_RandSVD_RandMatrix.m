D = 2e3;

A = randn(D, D)+randn(D, D)*1i;
A = expm(A+A')/norm(expm(A+A'));
Rank = 500;

tic
[U, S, V] = RandSVD(A,Rank, 'N_Oversamples', floor(Rank/4), 'N_Subspace_Iters', 0);
toc 

tic
[U_, S_, V_] = svd(A, 'econ'); 
  S_ = diag(S_); 
  S_ = S_(1:Rank);
  U_ = U_(:,1:Rank);
  V_ = V_(:,1:Rank);
toc

tic
Diff = norm(A - U*diag(S)*V')/norm(A)
toc

tic
Diff = norm(A - U_*diag(S_)*V_')/norm(A)
toc

figure();
semilogy(S_, 'r-o'); hold on
semilogy(S, 'b-x');
