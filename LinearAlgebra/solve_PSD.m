function [X,j] = solve_PSD(A,b,Minv)
%solve_PSD solves a linear system by means of the preconditioned steepest 
%          descent algorithm.
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
j = uint32(0);
w = 1;
alpha = 1;
while norm(r) > err && j < limit
    w = A*z;
    alpha = z'*r/(z'*w);
    X = X + z*alpha;
    r = r - w*alpha;
    z = Minv*r;
    j = j+1;
end
end
