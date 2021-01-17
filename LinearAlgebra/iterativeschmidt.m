function [Alist,slist,Vlist] = iterativeschmidt(psi,threshold)
%iterativeschmidt performs the iterative Schmidt decomposition on the
%wavefunction psi
%
% Arguments:
%     
%     psi (complex matrix): Matrix of wavefunction 
%     threshold (real number): Bond dimension to consider
%          
% Returns: Alist, slist, Vlist (complex lists): Tensor decompostion of the
%          wavefunction.
d1 = size(psi,1);
d2 = size(psi,2);
n = length(size(psi))-2;

psi = reshape(psi,d1*d2,d2^(n-1)*d1);

%%
[U,s,V] = svd(psi);
V = V';

% Truncate
chi = length(s(s>threshold));
s= s(1:chi,1:chi);
U= U(:,1:chi);
V= V(1:chi,:);

%%
A = U*s;
A = reshape(A,d1,d2,chi);

Alist{1} = A;
slist{1} = s;
Vlist{1} = V;

%%
for i = 1:n-2
%     red = reshape(V,d2*chi,d2^(n-i-1)*d1);
    [U,s,V] = svd(reshape(V,d2*chi,d2^(n-i-1)*d1));
    V = V';
    % Truncate
    chi2 = length(s(s>threshold));
    s= s(1:chi2,1:chi2);
    U= U(:,1:chi2);
    V= V(1:chi2,:);
%     A = U*s;
    A = reshape(U*s,chi,d2,chi2);
    Alist{i+1} = A;
    slist{i+1} = s;
    Vlist{i+1} = V;
    chi = chi2;
end
Vlist{end} = reshape(Vlist{end},chi,d2,d1);
end