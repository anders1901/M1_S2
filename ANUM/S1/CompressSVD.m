function J = CompressSVD(I, k)
    %R�alise une compression SVD de I avec k valeur singuli�re
    %On peut d�commenter la ligne suivante et supprimer l'argument k 
    %pour que l'utilisateur rentre k sous forme d'input
    %k = input('Enter a value between 0 and '+string(size(I,1))+' : ');
    [U,S,V] = svd(I,'econ');
    J = 0;
    for u = 1:k
        J = J + S(u,u)*U(:,u)*V(:,u)';
    end
    imagesc(J);
    colormap gray;
end
