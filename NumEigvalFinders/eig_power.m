function [vec,val] = eig_power(input_matrix)
%eig_power computes the closest eigenvector and eigenvalue of a given
% matrix.
% 
% Arguments:
%     
%     input_matrix (2D complex Hermitian matrix): matrix for the 
%     eigenvalue problem;
%          
% Returns: a right eigenvector and the corresponding eigenvalue of a
% matrix.

% Error
eabs = 1e-11;

% First estimated eigenvectors
n = size(input_matrix,1);
q = rand(n,1);
q = q/norm(q);
lambda_new = q'*input_matrix*q;
lambda_old = lambda_new - 10*eabs - 1;

while abs(lambda_new-lambda_old) > eabs
    lambda_old = lambda_new;
    q = (input_matrix)*q;
    q = q/norm(q);
    lambda_new = q'*input_matrix*q;
end

vec = q;
val = lambda_new;
end

