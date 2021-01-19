function result=mydft(input_array)
% mydft computes the discrete Fourier transform. It returns the exact same
%       result as the integrated fft function does for 1D arrays. 
% Note that this function might be illustrative for educational purposes.
%
% Arguments:
%
%       input_array (1D array complex): data to transform
%
% Returns: 1D array complex, transformed data of the same shape as
% an input array.

%% Info on array and discretised variables definition

qhorizontal = 0;
% Check for shape of array, we want to restrict the array to be a vector.
L = size(input_array);
if L(1) == 1
    input_array = input_array.';
    N = L(2);
    qhorizontal = 1;
elseif L(2) ==1
    N = L(1);
else
    disp('Error in input array shape')
    return
end

h = 2*pi/N; % define real space spacing

x = 0:h:(2*pi-h); % define real space grid, the 2pi factor seen in the 
                  % mathematical expression is alredy included in the x 
                  % variable
k = 0:1:N-1; % define recrp. space grid

if qhorizontal == 1 % to return the same shape of array as the input.
    result = (exp(-1i.*k'*x)*input_array).'; % computes dft (vectorized)
else
    result = (exp(-1i.*k'*x)*input_array); % computes dft (vectorized)
end


end