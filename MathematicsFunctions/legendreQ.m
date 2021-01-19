function [Q] = legendreQ(n,m,vargin)
%legendreQ obtains the syms expression for the Associated Legendre 
% function of the second kind of integer and real degrees n, m. Then, if a
% numerical vector is given in the argument, the numeric evaluation of the
% values in that vector is given at the output.
%   It's obtained via the Bonnet's recursion relation and Rodrigues'
%   formula
% Inputs: n, m: orders of the function in the sense Q_n^m
% Optional: x, vector to evaluate numerically the function
%
% REFERENCES: https://dlmf.nist.gov/14.9
%             https://dlmf.nist.gov/14.6
% Last update: September 26, 2020
% Juan Antonio Fernandez de la Garza -- juanfernandezdlg@gmail.com


% Short input chekup;
mNeg = 0;
if nargin > 3; error('Too many input argumments'); end
if mod(m,1) ~= 0 || mod(n,1) ~= 0; error('Order parameters are not integer'); end
if n+m<0; error('n+m >= 0 for this function to work'); end
if m < 0; m = -m; mNeg = 1; end

% Initialize variables and syms functions
syms s
syms Q(s)
syms Qn1(s);
syms Qn2(s);

Qn2 = 1/2 * log((1+s)/(1-s));
Qn1 = s*Qn2-1;

% Bonnet's recursion
if n==0
    Q = Qn2;
elseif n ==1
    Q = Qn1;
else
    for k = 2:n-1
        Q = (2*k-1)/k * s * Qn1 - (k-1)/k * Qn2;
        Qn2 = Qn1;
        Qn1 = Q; 
    end
    Q = (2*n-1)/n * s * Qn1 - (n-1)/n * Qn2;
end

%Rodrigues' formula
Q = (-1)^m * (1-s^2)^(m/2) * diff(Q,m);

% For m < 0
if mNeg == 1
    Q = (-1)^m * factorial(n-m)/factorial(n+m) *Q;
end

if nargin == 3
    if isnumeric(vargin)
        Qh = matlabFunction(Q);
        ax = feval(Qh,vargin);
        Q=ax;
    else
        error('Third argument is not a numerical vector')
    end
end


end

