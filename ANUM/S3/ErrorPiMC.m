function ErrorPiMC(N, a)
    % N : int, nombre max d'itérations
    % a : float, borne de l'interval [-a,a]
    if nargin<2
        % Dans l'exercice on fixe a = 1
        a =1;
    end
    R = 1:N;
    % On va stocker les valeurs successives de notre ?
    PI = []
    for n = R
        pi_ = PiMC(n, false, a);
        %err = abs(pi-pi_);
        PI = [PI pi_];
    end
    fprintf("Moyenne pour a = %d: %f \n", a, mean(PI));
    plot(R, PI); 
    hold on;
    plot([1 N],[pi pi]);
    title(['Loi forte des grands nombre sur [0,',num2str(a),']^2'])
    xlabel("Nombre max d'itérations"); 
    ylabel("Valeur de \pi approchée");
    axis([0 N 2.28 4]);
end