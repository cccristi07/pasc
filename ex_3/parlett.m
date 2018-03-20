function [ F ] = parlett( T, f )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    [m,n] = size(T);
    F = zeros(m);
    
    for i = 1:n
        F(i,i) = f(T(i,i));
    end
    
    for p = 1:n-1
        for i = 1:n-p
            
            j = i+p;
            s = T(i,j)*(F(j,j) - F(i,i));
            
            if p > 1
                for k = i+1:j-1
                    s = s + T(i,k)*F(k,j) - F(i,k)*T(k,j);
                    % s = T(i, i+1:j-1)*F(i+1:j-1,j) - F(i,i+1:j-1)*T(i+1:j-1, j)
                    % scris vectorizat...
                end
            end
            
            F(i,j) = s/(T(j,j)- T(i,i));

        end

    end
    
end

