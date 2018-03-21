function [ Ac, Bc, Cc ] = siso( num, den)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

   
    I = eye(length(den));
    
    % forma standard controlabila
    Ac = zeros(length(den)-1);
    Ac(1,:) = -den(2:end);
    for i = 2:length(Ac)
        Ac(i,i-1) = 1;
    end
    
    Bc = I(:,1);
    Cc = num;
    
    
    

end

