function [val,numrot] = eig_cj(input_matrix)
% eig_cj computes the eigenvalues of a matrix with the cyclic Jacobi
% algorithm
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


% find list of coordinates to cycle through
A = triu(input_matrix,1);
index = find(A);
[ind1,ind2] = ind2sub(size(A),index);

k = 1;

while off > eabs
    p = ind1(k);
    q = ind2(k);
    % Find c and s
    if input_matrix(p,q) ~= 0
        tau = (input_matrix(q,q) - input_matrix(p,p))/(2*input_matrix(p,q));
        if tau >= 0
            t = -tau + sqrt(1+tau^2);
        else
            t = -tau - sqrt(1+tau^2);
        end

        c = 1/sqrt(1+t^2);
        s = t*c;

        % Vector operations to make operation input_matrix*J
        input_matrix(:,[p,q]) = cat(2,input_matrix(:,p)*c-s*input_matrix(:,q),...
                                    s*input_matrix(:,p)+c*input_matrix(:,q));

        % Vector operations to make operation J'*(...)
        input_matrix([p,q],:) = cat(1,input_matrix(p,:)*c-s*input_matrix(q,:),...
                                    s*input_matrix(p,:)+c*input_matrix(q,:));
    
        off = norm(input_matrix.*(ones(n)-eye(n)));
        numrot = numrot + 1;

    end
    
    if k < length(ind1)
        k = k + 1;
    else
        k = 1;
    end

end

val = diag(input_matrix);

end

