x0 = [-1 -1].';
syms f2(x1,x2);
f2(x1, x2) = x1^2 + 2*x2^2;
[x1_, nb1, X1, E1] = Gradient(f2, x0, 0.1, 1e-16);
[x2_, nb2, X2, E2] = Gradient(f2, x0, 0.5, 1e-16);
[x3_, nb3, X3, E3] = Gradient(f2, x0, 1, 1e-16);
figure(1)
plot(0:nb1, E1); hold on;plot(0:nb2, E2);xlabel("Nombre d'itérations"); ylabel("Norme de l'erreur");
figure(2)
plot(0:nb3, E3);xlabel("Nombre d'itérations"); ylabel("Norme de l'erreur");

x0 = [-1 1.2].';
[T, t] = Wolfe(f, x0)