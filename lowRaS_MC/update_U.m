function [U,Hx] = update_U(x,V,Q,mu)
    Hx = hankel(x(1:floor(length(x)/2)),x(floor(length(x)/2):end));
%     R = rank(Hx);
    I = eye(size(V,2));
    U = mu * (Hx + 1/mu*Q ) * V / (I + mu * (V')*V);
end