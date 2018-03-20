function [ X ] = lyap_contR( A, C,varargin )
% Functie care calculeaza sol ec lyapunov continua 
% folosind exclusiv aritmetica reala
% A'X + XA - C  = 0 
% U * S' * U' * X + X * U * S * U' - C = 0
% S * U' * X * U + U' * X * U * S' - U' * C * U = 0
% S' * Y + Y * S + D = 0

%%
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
    
    if find( abs(eigA) <= 1e-10)
        disp('Avem valori proprii nule!!!')
        return
    end
    
    % obtinem forma Schur 
    [U, S] = schur(A, 'real');
       
    [m,~] = size(A);
    I = eye(m); 
    
    D = U' * C * U ; 
    
    
    Y = zeros(m);
    k = 1;
    
    
    if symmetricalC ~= true
        
        while ( k <= m )

            if k == m
                Y(:, m) = (S' + I * S(m,m)) \ ( D(:, m) - Y(:, 1:m-1) * S(1:m-1, m));
                break
            end

            if ( abs(S(k+1,k)) <= 1e-10)
                Y(:,k) = (S' + I*S(k,k)) \ ( D(:,k) - Y(:,1:k-1) * S(1:k-1,k));
                k=k+1;  
            else
                % avem bloc schur

                R = [S' + I * S(k,k), S(k+1,k) * I ;

                        S(k,k+1) * I, S'  + I * S(k+1,k+1)];

                d = [ D(:,k) - Y(:, 1:k-1) * S(1:k-1, k) ;
                    D(:,k+1) - Y(:, 1:k-1) * S(1:k-1, k+1)] ;

                y_aux = R \ d; 
                Y(:,k) = y_aux(1:m);
                Y(:,k+1) = y_aux(m+1:end);
                k=k+2;
            end
        end
        
    else
        while ( k <= m )
            if k == m
                 y = (S' + I * S(m,m)) \ ( D(:, m) - Y(:, 1:m-1) * S(1:m-1, m));
                 Y(m ,m) = y(end);
                 break
            end

            if ( abs(S(k+1,k)) <= 1e-10)
               % partitionam matricea R
              Rj = (S' + S(k, k) * eye(m));
              R21 = Rj(k:m, 1:k-1);
              R2 = Rj(k:m, k:m);
        
              f = D(k:m, k) - Y(k:m, 1:k-1) * S(1:k-1, k);
        
              yj2 = R2\(f - R21*Y(1:k-1, k));
        
              Y(k:m, k) = yj2;
        
              % simetrizam matricea Y
              Y(k, k:m) = yj2'; 
              k = k+1;
            else
            % avem bloc schur
            

             bM  = [S' + I * S(k,k), S(k+1,k) * I ;

                      S(k,k+1) * I, S'  + I * S(k+1,k+1)];
             
                  
             % aplicam factorizarea QR matricei bloc bM
             [Q, R] = qr(bM);
             

             d = [ D(:,k) - Y(:, 1:k-1) * S(1:k-1, k) ;
                   
                        D(:,k+1) - Y(:, 1:k-1) * S(1:k-1, k+1)] ;
             d =  Q' * d;
             
             
             
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
             size(y1);
             Y(k:m,k) = y1n;
             Y(k,k:m) = y1n;
             
             Y(k+1:m,k+1) = y2n(2:end);
             Y(k+1,k+1:m) = y2n(2:end);
             k=k+2;
            end
        end
            
        
        
    end
    
    
    X = U * Y * U';
    
    disp(norm( A'*X + X*A - C))



