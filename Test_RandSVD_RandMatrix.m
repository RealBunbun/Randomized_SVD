A = randn(1e3, 1e3);
A = A+A';
Rank = 100;

tic
[U, S, V] = RandSVD(A,Rank);
toc 

tic
[U_, S_, V_] = svd(A, 'econ'); 
  S_ = diag(S_); 
  S_ = S_(1:Rank);
  U_ = U_(:,1:Rank);
  V_ = V_(:,1:Rank);
toc

Diff = norm(A - U*diag(S)*V')/norm(A)
Diff = norm(A - U_*diag(S_)*V_')/norm(A)


figure();
semilogy(S_, 'r-o'); hold on
semilogy(S, 'bx');
