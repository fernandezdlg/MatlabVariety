function [val,numrot] = eig_j(input_matrix)
% eig_j computes the eigenvalues of a matrix with the classical Jacobi
% algor
% 
% Arguments:
% 
%     input_matrix (2D real symmetric matrix): matrix for the eigenvalue
%     problem;
% 
% Returns: val, an array with eigenvalues.
%          numrot, the number of rotations required.

% Initialization
n = size(input_matrix,1);
eabs = 1e-5*n^2;
off = norm(input_matrix.*(ones(n)-eye(n)));
numrot = 0;


while off > eabs
    % triu ensures p < q
    A = reshape(triu(input_matrix,1),n^2,1);
    [~,loc] = max(abs(A));
    [p,q] = ind2sub(size(input_matrix),loc);
    
    % I consider the if that was in the slides to be unnecessary as I pick
    % the max abs of the off-diag terms
    % Find c and s
    tau = (input_matrix(q,q) - input_matrix(p,p))/(2*input_matrix(p,q));
    if tau >= 0
        t = -tau + sqrt(1+tau^2);
    else
        t = -tau - sqrt(1+tau^2);
    end

    c = 1/sqrt(1+t^2);
    s = t*c;

%     % Matrix operations approach
%     J = eye(n);
%     J(p,p) = c;
%     J(q,q) = c;
%     J(p,q) = s;
%     J(q,p) = -s;
%     input_matrix = J'*input_matrix*J;
    % Vector operations to make operation input_matrix*J
    input_matrix(:,[p,q]) = cat(2,input_matrix(:,p)*c-s*input_matrix(:,q),...
                                s*input_matrix(:,p)+c*input_matrix(:,q));
  
    % Vector operations to make operation J'*(...)
    input_matrix([p,q],:) = cat(1,input_matrix(p,:)*c-s*input_matrix(q,:),...
                                s*input_matrix(p,:)+c*input_matrix(q,:));
    
    off = norm(input_matrix.*(ones(n)-eye(n)));
    numrot = numrot + 1;
end

val = diag(input_matrix);

end

