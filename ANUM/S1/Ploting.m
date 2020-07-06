function Ploting(a, b, eps, I)
    %Fonction qui permet d'afficher le taux de compression entre des bornes a
    %et b avec un pas eps sur une image I
    L = zeros(1, (b-a)/eps);
    i = 1;
    for r = a:eps:b
        J = CompressSVD(I, r);
        t = TauxCompress(I,J);
        L(i) = t;
        i = i+1;
    end
    plot(a:eps:b, log(L))
    title('Evolution du taux de compression selon k')
    xlabel('k : Nombre de valeurs singulières gardées')
    ylabel('\tau : log( Taux de compression )')
end