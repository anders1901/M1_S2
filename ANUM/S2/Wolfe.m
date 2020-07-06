function [T, t] = Wolfe(f, x0)
    t = 1;
    d = - gradient(f);
    %syms x1;
    %syms x2;
    syms x
    syms h(x);
    d_ = d(x0(1), x0(2));
    h(x) = f(x0(1)+x*d_(1),x0(2)+x*d_(2));
    H = diff(h);
    tg = 0;
    td = inf;
    
    m1 = 0.1;
    m2 = 0.9;
    while true
        if h(t)<=h(0) + m1*t*H(0)
            if H(t) >= m2*H(0)
                break
            else
                tg = t;
                if (td == inf) == 1
                    t = 10*tg; 
                else
                    t = (td + tg)/2;
                end
            end
        else 
            td = t;
            if (td == inf) == 1
               t = 10*tg; 
            else
                t = (td + tg)/2;
            end
        end
    end
    T = [tg td];
end