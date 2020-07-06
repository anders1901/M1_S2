function [x, nb, x_, E] = Gradient(f, x0, alpha, eps)
    err=10;
    df = gradient(f);
    xt = x0;
    x = x0;
    x_ = [x0];
    nb = 0;
    E = [0.0001];
    while(err>=eps & nb<=159)
        %fprintf("Etape : %d\n",nb)
        nb = nb+1;
        x = x - alpha*feval(df,x(1), x(2));
        err = norm(x-xt);
        E = [E err];
        xt = x;
        x_=[x_ x]; 
    end
end