function [ X ] = qr_lyap_contR( A, C )
% Functie care calculeaza sol ec lyapunov continua 
% folosind exclusiv aritmetica reala
% A'X + XA - C  = 0 (1)
% U * S' * U' * X + X * U * S * U' - C = 0
% S * U' * X * U + U' * X * U * S' - U' * C * U = 0
% S' * Y + Y * S + D = 0 (2)
% R' * Q' * Y + Y * Q * R - D = 0
% Q = Q' deci 
% R' * T * T * R - D = 0 (3)
% T = Y * Q
    
    % verificam matricea
    
    eigA = eig(A);
    
    if find( abs(eigA) <= 1e-10)
        disp('Avem valori proprii nule!!!')
        return
    end
    
    % obtinem forma Schur Complexa
    [U, S] = schur(A, 'real');
    
    [Q, R] = qr(S);
       
    [m,~] = size(A);
    I = eye(m); 
    
    D = U' * C * U ; 
    
    T = zeros(m);
    
    % aflam solutia eq (3)
    for i = 1:m
        
        T(:, i) = (R' + R(i, i)*I) \( D(:, i) - T(:, 1:i-1) * R(1:i-1, i));
        
    end
    
    issymmetric(C)
    % T = Y * Q = U' * X * U * Q
    Y = T * Q;
    X = U * Y * U';
    
    disp(norm( A' * X + X * A - C))

end

