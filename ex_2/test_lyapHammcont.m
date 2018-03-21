dim = 120;
A = rand(dim);
A = A*A';
A = A + dim*eye(dim);
A = -A;

C = rand(dim);
C = C*C';
Lc = chol(C);

Lx = lyap_contHamm(A,Lc);


