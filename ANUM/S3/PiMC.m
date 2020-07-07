function pi_ = PiMC(N, plot_, a)
    % N : int, nombre max d'it�rations
    % plot_ : logical, option pour afficher ou non rendu graphique
    % a : float, rayon 
    if nargin<3
        % Dans l'exercice on fixe de le quart de cercle unit�
        a = 1;
    end
    % On tire des nombre al�atoirs entre [0,a]
    x = a*rand(N,1);
    y = a*rand(N,1);
    
    % Formule du disque
    r = sqrt(x.^2+y.^2) ;
    % On obtient des logicals 
    inside = r<=a;
    outside = r>a;
    % plot 
    if plot_== true
        figure(1)
        plot(x(inside),y(inside),'b.');
        hold on
        plot(x(outside),y(outside),'r.');
        title(['Calcul approch� de \pi sur [0,',num2str(a),']^2'])
    end
    % Notre valeur de Pi
    pi_ = 4*sum(inside)/N;  
end