clear
clc
A = rand(10);


disp('test Taylor')
disp(norm(expTaylor(A,2) - expm(2*A)))

disp('test Pade')
disp(norm(expPade(A,2) - expm(2*A)))

% test 1/2
disp('Test 1/2')
disp(norm(scale_sq(A,2) - expm(2*A)))