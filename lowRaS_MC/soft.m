function a = soft(y, T)
% Soft-threshold function (for real or complex data)
% a = soft(y, T)
%
% INPUT
%    y : data
%    T : threshold
%
% If y and T are both multidimensional, then they must be of the same size.

a = max(1 - T./abs(y), 0) .* y;
