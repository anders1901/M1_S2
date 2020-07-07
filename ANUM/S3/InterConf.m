function [I1 I2] = InterConf(N, d)
    % N : int, nombre max d'itérations
    % d : int, dimension
    % Intervalle de confiance a 95%
    R = 1:N;
    % On va stocker les valeurs successives de notre Volume
    V = [];
    for n = R
        v = Volume(n, d);
        V = [V v];
    end
    mu = mean(V);
    sig = var(V);
    I1 = mu - 1.96*(sig/N)^0.5;
    I2 = mu + 1.96*(sig/N)^0.5;
    fprintf("Intervalle de confiance construit a partir de %d simulations.\n",N);
    fprintf("[%.6f, %.6f] \n",I1,I2);
end
