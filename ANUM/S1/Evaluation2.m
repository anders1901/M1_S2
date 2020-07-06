function res = Evaluation2(n, m, I)
    if n == m
        res = I;
    elseif n>m
        res = -1;
    else
        res = (I + exp(-1)*(n+m+1)) / (nchoosek(n+m,n)*factorial(m));
    end
 
end