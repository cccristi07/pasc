function [ F ] = scale_sq( A, t )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

    
    [n, ~] = size(A);
    m = 0;
    
    nA = norm(A);
    
    while t*nA >= 1
        t = t/2;
        m = m+1;
    end
    
    F = expm(t*A);
    while m >= 1
        F = F^2;
        m = m-1;
    end
    
end

