x0 = [-1 -2].';

figure(1)
[x, nb, x_, E] = Newton_Rosenbrock(f, x0, 1e-16);

figure(2)
plot(0:nb, E);ylabel("Norme de l'erreur");
