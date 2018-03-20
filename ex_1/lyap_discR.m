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
    [U, S] = schur(A, 'complex');
       
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
        disp('Eroare?')
        while ( k <= m )

            if k == m 
                d = D(k, k) - S(k, k)' * Y(k, 1:k-1) * S(1:k-1, k);
                Y(k, k) = (S(k,k)' * S(k, k) + 1) \ d;
                break
            end

            if ( abs(S(k+1,k)) <= 1e-10)
                   
                d = D(k:m, k) - S(k:m, k:m)' * Y(k:m, 1:k-1) * S(1:k-1, k);
                R = (S' * S(k,k) + eye(m));
                R21 = R(k:m, 1:k-1);
                R2 = R(k:m, k:m);
                
                y = R2 \ ( d - R21 * Y(1:k-1, k));
                % simetrizam pe Y
                Y(k:m, k) = y;
                Y(k, k+1:m) = y(2:end)';
                k = k + 1;
            else
                
                % avem bloc schur
                A_ = S' * S(k,k) + eye(m);
                B_ = S' * S(k, k+1);
                C_ = S' * S(k+1, k);
                D_ = S' * S(k+1, k+1) + eye(m);
                
                A21 = A_(k:m, 1:k-1);
                A2 = A_(k:m, k:m);
                
                B21 = B_(k:m, 1:k);
                B2 = B_(k:m, k+1:m);
                
                C21 = C_(k:m, 1:k-1);
                C2 = C_(k:m, k:m);
                
                D21 = D_(k:m, 1:k);
                D2 = D_(k:m, k+1:m);
                
                R = [ A21 A2 B21 B2; C21 C2 D21 D2];
               

                d = [ D(k:m,k) - Y(k:m, 1:k-1) * S(1:k-1, k) ;
                    D(k:m,k+1) - Y(k:m, 1:k-1) * S(1:k-1, k+1)] ;

                y_aux = R \ d; 
                Y(k:m,k) = y_aux(k:m);
                Y(k, k+1:m) = y_aux(k+1:m)';
                
                Y(k+1:m,k+1) = y_aux(m+k+1:end);
                Y(k, k:m) = y_aux(m-k+2:end)';
                k=k+2;

            end
        end

        
        
    end
    
    X = U * Y * U';

    disp(norm( A' * X * A + X - C))

end

