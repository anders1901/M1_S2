function res = Evaluation(n)
    if n == 0
        res = 1 - exp(-1);

    else
        res = -exp(-1) + n*Evaluation(n-1);
    end
end
    
    
