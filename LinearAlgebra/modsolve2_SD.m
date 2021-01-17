function [relerr,enorm] = modsolve2_SD(A,b,N)
%modsolve2_SD solves a linear system by means of the steepest descent 
%            algorithm and outputs every relative error and error norm
% 
% Arguments:
%     
%     A (complex matrix): Matrix of the linear system 
%     b (complex vector): Vector of the linear system
%     N (integer number): Number of iterations to be performed
%          
% Returns: relerr (real vector): vector with the relative error at each
%               iteration.
%          enorm (real vector): vector with the ||e_i||_A at each
%               iteration.

X = rand(length(A),1);
relerr = zeros(N,1);
enorm = zeros(N,1);
X0 = A\b;
r = b-A*X;
relerr(1) = norm(X0 - X);
enorm(1) = sqrt((X-X0)'*A*(X-X0));

for k=2:N
    X = X + r'*r./(r'*A*r)*r;
    r = b-A*X;
    relerr(k) = norm(X0 - X);
    enorm(k) = sqrt((X-X0)'*A*(X-X0));
end
relerr = relerr./norm(X0);
end
