%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vandermonde                                          %
% Edit by Min Ding 30/08/2019                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = vander_add (v,N)
    if size(v,1)~= 1
        v = v';
    else
        v = v;
    end
    M = length(v);
    A = repmat(v,N,1);
    A(1,:) = 1;
    A = cumprod(A);
end