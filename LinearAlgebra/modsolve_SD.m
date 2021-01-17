function [X] = modsolve_SD(A,b,N)
%modsolve_SD solves a linear system by means of the steepest descent 
%            algorithm and outputs every estimation of X.
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
err = 1e-14;
err = err*norm(b);
r = b-A*X(:,1);
for k=2:N
    X(:,k) = X(:,k-1) + r'*r./(r'*A*r)*r;
    r = b-A*X(:,k);
end
end
