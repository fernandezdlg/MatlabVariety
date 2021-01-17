function [X,j] = solve_CG(A,b)
%solve_CG solves a linear system by means of the conjugate gradient method.
%
% Arguments:
%     
%     A (complex matrix): Matrix of the linear system 
%     b (complex vector): Vector of the linear system
%          
% Returns: X (complex vector): Solution of the linear system
%          j (real number): Number of iterations performed
limit = uint32(10000000);
err = 1e-14;
err = err*norm(b);
X = rand(length(A),1);
r = b-A*X;
d = r;
alpha = 1;
ri = r;
j = uint32(0);
while norm(r) > err  && j < limit
    ri = r;
    alpha = r'*r/(d'*A*d);
    X = X + alpha*d;
    r = r - alpha*A*d;
    d = r + r'*r/(ri'*ri)*d;
    j = j + 1;
end
end
