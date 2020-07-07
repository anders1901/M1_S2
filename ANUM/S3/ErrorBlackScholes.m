function ErrorBlackScholes(N, K, param)
    % N : int, nombre max d'itérations
    % K : float, prix 
    % param : array(float,2), mu et sigma de N(mu,sigma)
    if nargin <3
        if nargin<2
            %Dans l'exercice on suppose K = 1
            K = 1;
        end
        %On travail sur une loi normale centrée réduite
        param = [0 1];
    end
    R = 1:N;
    % On va stocker les valeurs successives de notre ?
    C_ = [];
    P_ = [];
    c_ = exp(1/2)*normcdf(1)- 0.5;
    p_ = 0.5 - exp(1/2)*normcdf(-1);
    for n = R
        [c p] = BlackScholes(n);
        %err = abs(pi-pi_);
        C_ = [C_ c];
        P_ = [P_ p];
    end
    fprintf("Moyenne pour C : %f \n", mean(C_));
    fprintf("Moyenne pour P : %f \n", mean(P_));
    %figure(1)
    plot(R, C_); 
    hold on;
    plot(R, P_); 
    hold on;
    plot([1 N],[c_ c_]);
    hold on;
    plot([1 N],[p_ p_]);
    title(['Loi forte des grands nombre pour N = ',num2str(N),'.'])
    xlabel("Nombre max d'itérations"); 
    ylabel("Valeur de C et P approchée");
    axis([0 N 0 1.6]);
end