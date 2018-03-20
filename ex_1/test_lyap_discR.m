A = rand(4);
C = rand(4);

C = C +  C';

%lyap_discR(A, C);

lyap_discR(A,C,false);