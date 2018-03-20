function [ X ] = lyap_discR( A, C, varargin )
% Functie care calculeaza sol ec lyapunov continua 
% folosind exclusiv aritmetica reala
% A'X A - X = C 
% U * S' * U' * X * U * S * U' + X = C 
% S' * U' * X * U * S * U' + U' * X * U = U' * C * U
% S' * Y * S + Y  = D 
% S
%
%
%% Verificam daca matricea e simetrica sau nu

if nargin == 3
    symmetricalC = varargin{1};
elseif nargin == 2
    symmetricalC = false;
elseif nargin > 3
    error('Prea multe argumente!')
end
%%

    
    % verificam matricea
    
    eigA = eig(A);
    
     if intersect(eigA, -1./eigA)
        disp('Matricea A nu este valida')
        return
    end
    
    % obtinem forma Schur Complexa
    [U, S] = schur(A, 'real');
       
    [m,~] = size(A);
    I = eye(m); 
    
    D = U' * C * U ; 
    
    
    Y = zeros(m);
    k = 1;
    
    if symmetricalC ~= true
        
        while ( k <= m )
            
            if k == m 
                
                d = D(:, k) - S' * Y(:, 1:k-1) * S(1:k-1, k);
                Y(:, k) = (S' * S(k, k) + eye(m)) \ d;
                break


            elseif ( abs(S(k+1,k)) <= 1e-10)
                
                d = D(:, k) - S' * Y(:, 1:k-1) * S(1:k-1, k);
                Y(:, k) = (S' * S(k, k) + eye(m)) \ d;
                k = k +1;
                
            elseif ( abs(S(k+1,k)) >=  1e-10)
                
                d1 = D(:, k) - S' * Y(:, 1:k-1) * S(1:k-1, k);
                d2 = D(:, k+1) - S' * Y(:, 1:k-1) * S(1:k-1, k+1);
                dk = [ d1; d2];
                R = [ S' * S(k,k) + I, S' * S(k+1, k);
                    S' * S(k, k+1), S' * S(k+1, k+1) + I];

                y = R \ dk;

                Y(:, k) = y(1:m);
                Y(:, k+1) = y(m+1:end);
                k=k+2;  


            end
        end


    else
        
        while ( k <= m )

            if k == m 
                d = D(:, k) - S' * Y(:, 1:k-1) * S(1:k-1, k);
                y = (S' * S(k, k) + eye(m)) \ d;
                Y(m,m) = y(end);
                
                break
            end

            if ( abs(S(k+1,k)) <= 1e-10)
                   
                d = D(:, k) - S' * Y(:, 1:k-1) * S(1:k-1, k);
                R = (S' * S(k,k) + eye(m));
                R21 = R(k:m, 1:k-1);
                R2 = R(k:m, k:m);
                d = d(k:m);
                y = R2 \ ( d - R21 * Y(1:k-1, k));
                % simetrizam pe Y
                Y(k:m, k) = y;
                Y(k, k+1:m) = y(2:end)';
                k = k + 1;
            else
                % avem bloc schur
                A_ = S' * S(k,k) + eye(m);
                B_ = S' * S(k+1, k);
                C_ = S' * S(k, k+1);
                D_ = S' * S(k+1, k+1) + eye(m);
                
                bMat = [A_ B_; C_ D_];
                
                [Q, R] = qr(bMat);
                
                d = Q' * [ D(:,k) - S' * Y(:, 1:k-1) * S(1:k-1, k) ;
                    D(:,k+1) - S' * Y(:, 1:k-1) * S(1:k-1, k+1)] ;
                
                d1 = d(1:m);
                d2 = d(m+1:end);
             
                d2n = d2(k:m);
                d1n = d1(k:m);
             
                d2c = d2(1:k-1);
                d1c = d1(1:k-1);
             
             
                len_nec = length(d2n);
                R2 = R(end-len_nec+1:end, end-len_nec+1:end);
             
                y2c = Y(1:k-1,k+1);
                y2n = R2 \ d2n;
             
             
                R1 = R(k:m, k:m);
                R12 = R(k:m, m+1:m+k-1);
                R13 = R(k:m, m+k:end);
             
                y1n = R1 \ (d1n - R12*y2c - R13*y2n);
             
                y_aux = R \ d;
                y1 = y_aux(1:m);
             
                %y2 = y_aux(m+1:end);
      
                Y(k:m,k) = y1n;
                Y(k,k:m) = y1n;
             
                Y(k+1:m,k+1) = y2n(2:end);
                Y(k+1,k+1:m) = y2n(2:end);
%                 y_aux = R \ d; 
%                 Y(k:m,k) = y_aux(k:m);
%                 Y(k, k+1:m) = y_aux(k+1:m)';
%                 
%                 Y(k+1:m,k+1) = y_aux(m+k+1:end);
%                 Y(k, k:m) = y_aux(m-k+2:end)';
                k=k+2;

            end
        end

        
        
    end
    
    X = U * Y * U';

    disp(norm( A' * X * A + X - C))

end

