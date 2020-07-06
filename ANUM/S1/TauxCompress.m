function t = TauxCompress(I,J)
    %Calcule le taux de compression de 
    % - J l'image compressée
    %par rapport à 
    % - I l'image originelle 
    [U,S,V] = svd(I,'econ');
    [u,s,v] = svd(J,'econ');
    t = norm(s, "fro")/norm(S, "fro");
end 