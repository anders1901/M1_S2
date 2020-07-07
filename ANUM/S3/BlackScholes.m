function [C P] = BlackScholes(N, K, param)
    % N : int, nombre max d'itérations
    % K : float, prix 
    % param : array(float,2), mu et sigma de N(mu,sigma)
    if nargin <3
        if nargin<2
            %Dans l'exercice on suppose K = 1
            K = 1;
        end
        %On travail sur une loi normale centrée réduite
        param = [0 1];
    end
    % On obtient des nombres tirés selon N(mu,sigma)
    N_ = normrnd(param(1), param(2), N, 1);
    %r = exp(N_) - 1;
    % On obtient des logicals 
    inside = N_ >= 0;
    outside = N_ < 0;
    C = sum(exp(N_(inside)) - 1)/N;
    P = sum(1 - exp(N_(outside)))/N;
end