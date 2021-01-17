function result = mydftinverse(input_array)
% mydftinverse computes the discrete inverse Fourier transform.
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
% Check for shape of array
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
    result = (h/(2*pi)).*(exp(1i.*k'*x)*input_array).'; % computes idft
else
    result = (h/(2*pi)).*(exp(1i.*k'*x)*input_array); % computes idft
end


end