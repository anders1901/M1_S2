function [x, nb, x_, E] = Newton_Rosenbrock(f, x0,eps) 
    err=1;
    df = gradient(f);
    H = hessian(f);
    xt = x0;
    x = x0;
    x_ = [x0];
    nb=0;
    E = [0.0001];
    while (nb <5) 
        nb = nb+1;
        %x = x - feval(df,x(1), x(2))feval(f,x(1), x(2)); 
        x = x - inv(feval(H,x(1), x(2)))*feval(df,x(1), x(2)); 
        err=norm(x-[1;1]);
        E = [E err];
        xt = x;
        x_=[x_ x];
    end %while
    scatter(x_(1,:), x_(2,:)); 
    hold on;
    ezcontour(f, [-1.5;2;-3;3])
end