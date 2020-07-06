syms x1;
syms x2;
syms f(x1,x2);
f(x1,x2) = (x1 - 1)^2 + 100*(- x1^2 + x2)^2;
fprintf("f(x1,x2) = %s \n",f)
g = jacobian(f,[x1,x2]);
H = hessian(f, [x1,x2]);

%fprintf("Le gradient de f est : g(x1,x2) = %s \n",g)
disp(" ")
disp("- Le gradient de f est : ")
disp("g(x1, x2) =")
disp(g(x1,x2).')

%fprintf("La matrice Hessienne de f est : H(x1,x2) = %s \n",H)
disp("- La matrice Hessienne de f est : ")
disp("H(x1, x2) =")
disp(H(x1,x2))

x_ = solve(g(x1,x2)==0);
%s = sprintf('Le seul point critique de f est : x_ = (x1_, x2_) =  (%d, %d) ', x_.x1, x_.x2);
%disp(s)
fprintf("- Le seul point critique de f est : \nx_ = (x1_, x2_) =  (%d, %d)\n\n", x_.x1, x_.x2)

s = sprintf('Or det(H(1, 1)) = %d > 0', det(H(1,1)));
disp(s)

disp("Donc c'est bien un extremum local.")

s = sprintf('De plus tr(H(1, 1)) = %d > 0', trace(H(1,1)));
disp(s)

disp("Donc c'est bien un minimum local.")
