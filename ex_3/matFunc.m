function [ F ] = matFunc( A, f )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    
    [U, S] = schur(A, 'complex');
    
    T = parlett(S, f);
    
    F = U*T*U';
    
    F = real(F);

end

