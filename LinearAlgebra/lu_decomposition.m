function [L,U,P] = lu_decomposition(A)

% lu_decomposition computes the LU decomposition with pivoting
% for a generally complex square matrix A
% Return lower L, upper U and permutation P matrices such that P*A = L*U

n = issquare(A);

L = zeros(n);
P = eye(n);

for k = 1:n-1
    P1 = eye(n);
    [~,c] = max(abs(A(k:n,k))); % largest element index
    if c ~= k                % for robustness, diagonal term is abs larger
        A([k,c+k-1],:) = A([c+k-1,k],:);
        P([k,c+k-1],:) = P([c+k-1,k],:);
        L([k,c+k-1],:) = L([c+k-1,k],:);        
    end
    L(k:n,k) = A(k:n,k)./A(k,k);
    P1(:,k) = 2*P1(:,k)-L(:,k);
    A = P1*A;
end

L(end,end)=1;
U = A;

end