function V = Volume(N, d, plot_, a)
    % N : int, nombre max d'itérations
    % d : int, dimension
    % plot_ : logical, option pour affihcer
    % a : float, borne de l'interval [-a,a]
    if nargin<4
        if nargin < 3
            % On n'affiche pas tout le temps 
            plot_ = false;
        end
        % Dans l'exercice on fixe a = 1
        a = 1;
    end
    % On récupére des nombre ~ U([-a,a])
    X = 2*a*rand(N,d) - a;
    % X2 ser pour r (partie somme de la norme 2
    X2 = X.^2;
    %norme 
    r = sqrt(sum(X2(:,:),2));
    % On obtient des logicals 
    inside = r<=a;
    outside = r>a;
    % Notre valeur du volume 
    V = (2*a)^d*sum(inside)/N;
    if plot_ == true
        if d==3
            % Si d = 3 on affiche le résultat
            figure(1)
            plot3(X(inside,1),X(inside,2),X(inside,3),'b.');
            hold on
            plot3(X(outside,1),X(outside,2),X(outside,3),'r.');
        end
        if d == 2
            % Si d = 2 on affiche le résultat
            figure(1)
            plot(X(inside,1),X(inside,2),'b.');
            hold on
            plot(X(outside,1),X(outside,2),'r_');
        end
    end
end