function [x, X, E, iter] = GradientConjugue(A, y, x0, maxit, eps)
x = x0;
xt = x;
X = [];

err = 1;
E = [err];

n = size(A,1);
r = y - A*x;
p = r;
xs = A\y;
iter = 1;
while (iter <= maxit && err>eps)
  d = dot(A*p, p);
  %alpha = dot(r,p)/d;
  alpha = dot(r,p)/d;
  xt = x + alpha*p;
  r = r - alpha*A*p;
  
  beta = - dot(A*r, p)/d;
  p = r + beta*p;
  
  err = norm(x - xt);
  E = [E; err];
  x = xt;
  X = [X; x];
  iter = iter + 1;
end
