function [X,j] = solve_PCG(A,b,Minv)
%solve_PCG solves a linear system by means of the preconditioned conjugate
%          gradient method.
% 
% Arguments:
%     
%     A (complex matrix): Matrix of the linear system 
%     b (complex vector): Vector of the linear system
%     Minv (complex matrix): Inverse of the preconditioning matrix
%
% Returns: X (complex vector): Solution of the linear system
%          j (real number): Number of iterations performed
limit = uint32(10000000);
err = 1e-14;
err = err*norm(b);

X = rand(length(A),1);
r = b-A*X;
z = Minv*r;
p = z;

alpha = 1;
ri = r;
zi = z;
j = uint32(0);
while norm(r) > err  && j < limit
    alpha = ri'*zi/(p'*A*p);
    X = X + alpha*p;
    r = r - alpha*A*p;
    z = Minv*r;
    p = z + r'*z/(ri'*zi)*p;
    j = j + 1;
    ri = r;
    zi = z;
end
end
