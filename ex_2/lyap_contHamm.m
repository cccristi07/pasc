function [ Lx ] = lyap_contHamm( A, C )
% A'X + XA  + C = 0
% U*S'*U'*X + X*U*S*U' + C = 0
% S'* U'*X*U + U'*X*U*S + U'*C*U = 0
% S e bloc sup triungiulara
% 

    [m, n] = size(A);
    %desc schur reala a lui A
    [U, S] = schur(A, 'complex');
    
    % fact cholesky a lui C
    Lc = chol(U'*C*U);
    
    % prim pas
    j = 1;
    Lx = zeros(m);
    B = Lc;
    while j <= m
            % calculam factorul Q al matricii B(:, j)
            
            [Q, ~] = qr(B(:,j));
            
            %Q * R = Lc 
            
            % aplicam fact Q asupra lui B
            B(:, j:m) = Q' * B(:,j:m);
            
            Lx(j, j) = B(1, j) / sqrt(-2*S(j,j));
            
            fj = B(1,j)/Lx(j,j);
            
            r = -B(1, j+1:n)'*fj - S(j,j+1:n)'*Lx(j,j);
            u = (S(j+1:n,j+1:n)' + S(j,j)*eye(m-j))\r;
       
            Lx(j, j+1:end) = u';
            B(1,j+1:m) = B(1,j+1:m) - fj*Lx(j,j+1:m);
            
            j = j+1;
      
    end
    
    Lx = U * Lx * U';
    X = Lx' * Lx;
    norm(A'*X + X*A + C)
  
    
end

