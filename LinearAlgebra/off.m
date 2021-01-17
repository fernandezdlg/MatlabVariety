function res = off(A)
% off computes the norm of the off-diagonal elements of matrix A.
% 
% Arguments:
%     
%     A: 2D arbitrary matrix.
% 
% Returns: off-diagonal norm value.
n = size(A,1);
A = A - diag(diag(A));
A = reshape(A,n^2,1);

res = norm(A);


end

