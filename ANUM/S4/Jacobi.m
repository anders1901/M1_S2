function [x, X, E, iter] = Jacobi(A, y, x0, maxit, eps)
x = x0;
xt = x;
X = [];

D = diag(diag(A));
L = -tril(A,-1);
U = -triu(A,1);
%normf = norm(y);
err = 1;
E = [err];
xs = A\y;
iter = 1;
while (iter <= maxit && err>eps)
  % xt = x - D\(A*x + y);
  xt = D\(y+(L+U)*x);
  err = norm(x - xt);
  E = [E; err];
  x = xt;
  X = [X; x];
  iter = iter + 1;
end

