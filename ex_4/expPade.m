function [F ] = expPade( A, t )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
    
    [m, n] = size(A);
    
    eps = 1e-5;
    
    X = eye(m);
    N = eye(m);
    D = eye(m);
    a = t*norm(A);
    
    for p = 1:100
        
        if 8 * (factorial(p)^2) * ( t*norm(A))^(2*p+1)/(factorial(2*p) *factorial(2*p+1)) <= eps;
            break;
        end
    end
    
    for k = 1:p
        X = (p-k+1)/(k*(2*p-k+1)) *t*A*X;
        N = N + X;
        D = D + (-1)^k * X;
    end
    F = D\N;

end

