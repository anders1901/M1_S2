%format longE;
function res = Higham(x)
    for i=1:52
        x=sqrt(x);
    end
    for i=1:52
        x=x.^2;
    end 
    res = x;
end

