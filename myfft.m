function result = myfft(input_array)
% myfft computes the discrete Fourier transform with the radix-2 algorithm
% Note that this function might be illustrative for educational purposes.
%
% Arguments:
%
%   input_array (1D complex array, size 2^N for N>0): data to transform;
%
% Returns: 1D complex array, transformed data of the same shape as the
% input array.

%% Check for shape of array
horizontal = 0;
L = size(input_array);
if L(1) == 1
    input_array = input_array.';
    L = L(2);
    horizontal = 1;
elseif L(2) == 1
    L = L(1);
else
    error('Error in input array shape (not 1D)')
end

%% Check for 2^N length
Factors = factor(L);
if isempty(find(Factors>2,1))
    N = length(Factors);
else
    error('Input array not 2^N')
end

%% Bit reverse and definitions
input_array = bitrevorder(input_array);
h = 2*pi/L;
Ph = 0:h:(pi-h);
W = exp(-1j*Ph).';

%% Propagation
delta = 2;  % separation between crosses groups convenient variable

for layer = 1:N             % layer counter
    for m = 1:(2^(N-layer))   % crosses group counter
        %       Reminder of how to call G or H on a given cross:  
        %       input_array(((m-1)*delta+1):(m*delta -delta/2)))  % G
        %       input_array((m*delta -delta/2 + 1):(m*delta))) %H

        % Here H is multiplied by the Twiddle factors
        input_array((m*delta -delta/2 + 1):(m*delta)) = ...
            input_array((m*delta -delta/2 + 1):(m*delta)).* ...
            W(1:(2^(N-layer)):end); %  H turned to H'
        
        % Here G and H are combined to form the next layer array
        input_array(((m-1)*delta+1):(m*delta)) = ...
            cat(1,input_array(((m-1)*delta+1):(m*delta -delta/2)), ...
                input_array(((m-1)*delta+1):(m*delta -delta/2))) + ...
            cat(1,input_array((m*delta -delta/2 + 1):(m*delta)), ...
                -1.*input_array((m*delta -delta/2 + 1):(m*delta))); 
    end
    delta = delta*2;
end

if horizontal == 1
    result = input_array.';
else
    result = input_array;
end

end
