function F = Defloute(A, B, G, p)
    % A premiere partie de K
    % B deuxieme partie de K
    % G image floutée 
    % p nombre de valeurs singuliéres
    n = size(G, 2);
    
    [Ua,Sa,Va] = svd(A,'econ');
    [Ub,Sb,Vb] = svd(B, 'econ');
    Ghat = Ub.'*G*Ua;
    S = diag(Sb)*(diag(Sa)).';
    
    s = reshape(S,n^2,1);
    [ss,prm] = sort(s);
    prm = flipud(prm);
    ss = flipud(ss);
    ls = length(s);
    iprm = zeros(ls,1);
    iprm(prm) = [1:ls]';

    ssnew = [ss(1:p); zeros(length(ss)-p,1)];
    Snew = reshape(ssnew(iprm),n,n); 
    
    Fhat = Ghat./ S;
    Fnew = Fhat .* (Snew>0);
    F = Vb*Fnew*Va';
    imshow(F);
end