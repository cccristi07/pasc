A = rand(5);


% radical 
F = matFunc(A, @(t) sqrt(t));

disp( F*F - A);

% 1/z
F = matFunc(A, @(t) 1/t);
disp(A*F);

% ln(z)
F = matFunc(A, @(t) log(t));
disp(expm(F) - A);