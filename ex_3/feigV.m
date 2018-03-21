function [F] = feigV(A, f)

[eigVec, eigVal] = eig(A);

F = eigVec * f(eigVal) * eigVec^(-1);


end