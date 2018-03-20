dim = 5;

eigv = rand(dim, 1);

A = full(sprandsym(dim, 0.7, eigv));
C = rand(dim);
C = C*C';

lyap_discHamm(A, C);