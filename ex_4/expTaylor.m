function [ F ] = expTaylor( A, t )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
    [m, n] = size(A);
    tol = 1e-5;
    X = eye(m);
    F = eye(m);
    
    a = t*norm(A);
    
    for p = 1:100
        
        if a^(p+1) / ( factorial(p+1)) <= tol
            break
        end
    end
    
    for k = 1:p
        X = (1/k)*t*A*X;
        F = F + X;
    end
end

