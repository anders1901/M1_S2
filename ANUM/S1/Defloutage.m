function F = Defloutage(A, B, G, p)
    [Ua,Sa,Va] = svd(A,'econ');
    [Ub,Sb,Vb] = svd(B, 'econ');
    Ghat = Ub.'*G*Ua;
    S = diag(Sb)*(diag(Sa)).';
    Fhat = Ghat./ S.*S+lam^2);
    F = Vb*Fhat*Va.';
    
    U = kron(Ua,Ub);
    S = kron(Sa, Sb);
    V = kron(Va, Vb);
    g = reshape(G.', [size(G,1),1]);
    T = 0;
    for i = 1:p
        T = T+ U(:,i).'*V(:,i)*g(1:p)/S(i,i);
    end
    
    %F = (Vb/Sb)*Ub.'*G*(Ua/Sa)*Va.';
    %for i = 1:p
     
end