function [X,j] = solve_SD(A,b)
%solve_SD solves a linear system by means of the steepest descent algorithm.
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
j = uint32(0);
while norm(r) > err && j < limit
    X = X + r'*r./(r'*A*r)*r;
    r = b-A*X;
    j = j+1;
end
end
