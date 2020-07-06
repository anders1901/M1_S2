function l_ = PowerM(A, m, n)
    %Méthode de la puissance sur la matrice A 
    % Seuil max : m
    %dimension de x : n 
    x = rand(n,1);
    l=0;
    l_ = 100;
    %for i = 1:m
    while abs(l-l_)>=m
        y_ = A*x;
        x_ = y_./norm(y_);
        l=l_;
        l_ = x_.'*A*x_;
        %i = i+1
    end
end