function [x, X, E, iter] = GausseSeidel(A, y, x0, maxit, eps)
x = x0;
xt = x;
X = [];

err = 1;
E = [err];

n = size(A,1);
D = diag(diag(A));
L = -tril(A,-1);
U = -triu(A,1);
xs = A\y;
iter = 1;
while (iter <= maxit && err>eps)
  xt = (D-L)\(y + U*x);
  err = norm(x - xt);
  E = [E; err];
  x = xt;
  X = [X; x];
  iter = iter + 1;
end
