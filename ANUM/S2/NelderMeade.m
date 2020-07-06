function [X_, n] = NelderMeade(X0, f)
    X_ = X0;
    max_eval = 250;
    n = 0;
    while n<max_eval
        n = n+1;
        F_ = f(X_(:,1),X_(:,2));
        [B I] = sort(F_);
        N_ = size(F_,1);
        X_temp = [];
        for i = 1:N_
            X_temp(i,:) = X_(I(i),:);
        end
        X_ = X_temp;
        x_0 = mean(X_(1:end-1,:));
        x_N = X_(end-1,:);
        x_r = x_0 + (x_0 - X_(end,:));

        if f(x_r(1), x_r(2)) < f(x_N(1), x_N(2))
            x_e = x_0 + 2*(x_0 - X_(end,:));
            if f(x_e(1), x_e(2)) < f(x_r(1), x_r(2))
                X_(end,:) = x_e;
            else    
                X_(end,:) = x_r;
            end
        end
        if f(x_r(1), x_r(2)) > f(x_N(1), x_N(2))
            x_c = X_(end,:) + (x_0 - X_(end,:))/2;
            if f(x_c(1), x_c(2)) < f(x_N(1), x_N(2))
                X_(end,:) = x_c;
            else
                for i =2:N_
                    X_(i,:) = X_(1,:) + (X_(i,:) - X_(1,:))/2;
                end
            end
        end
    end
end