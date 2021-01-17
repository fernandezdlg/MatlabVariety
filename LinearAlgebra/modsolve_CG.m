function [X] = modsolve_CG(A,b,N)
%modsolve_CG solves a linear system by means of the conjugate 
%            gradient method and outputs every estimation of X.
%
% Arguments:
%     
%     A (complex matrix): Matrix of the linear system 
%     b (complex vector): Vector of the linear system
%     N (integer number): Number of iterations to be performed
%          
% Returns: X (complex vector): All estimates of the linear system
X = zeros(length(A),N);
X(:,1) = rand(length(A),1);
r = b-A*X(:,1);
d = r;
alpha = 1;
ri = r;
for k=2:N
    ri = r;
    alpha = r'*r/(d'*A*d);
    X(:,k) = X(:,k-1) + alpha*d;
    r = r - alpha*A*d;
    d = r + r'*r/(ri'*ri)*d;
end
end
