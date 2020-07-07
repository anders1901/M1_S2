function Erreur(eps, N, opt)

%Comparaison de l'erreur des méthodes 
%Jacobi
%Gauss-Seidel
%SOR
%Gradient conjugué 
%opt = 1 matrice poisson
%opt = -1 matrice binomial
%opt = 0 matrice binomial
Nmax = 1000; % nombre d'itérations max
m = 1:1:Nmax;
if opt==1
    A = gallery('poisson', N); 
    y = zeros(N*N, 1); 
    x0 = ones(N*N, 1); 
end
if opt==-1
    A = gallery('binomial', N); 
    y = randn(N, 1); 
    x0 = ones(N, 1); 
end
if opt==0
    A = randn(N,N);
    y = zeros(N, 1); 
    x0 = ones(N, 1); 
end
[x_jac, X_jac, E_jac, iter_jac] = Jacobi(A, y, x0, Nmax, eps);
[x_gs, X_gs, E_gs, iter_gs] = GausseSeidel(A, y, x0, Nmax, eps);
[x_sor, X_sor, E_sor, iter_sor] = Sor(A, y, x0, 1.8, Nmax, eps);
[x_gc, X_gc, E_gc, iter_gc] = GradientConjugue(A, y, x0, Nmax, eps);

clf;
semilogy(m(1:iter_jac), E_jac, 'r-', m(1:iter_gs), E_gs, 'b-', m(1:iter_sor), E_sor, 'c-', m(1:iter_gc), E_gc, 'g-', 'linewidth', 2);
disp(E_jac(1));
ylim([eps max([max(E_jac) max(E_gs) max(E_sor) max(E_gc)])]);
set(gca, 'fontsize', 18);
legend('Jacobi', 'GS', 'SOR', 'CG');
c = sprintf("Convergence pour eps = %.1e", eps);
title(c);
ylabel('Erreur');
xlabel("Nombre d'itérations");
