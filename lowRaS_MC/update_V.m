function V = update_V(Hx,U,Q,mu)
    I = eye(size(U,2));
    V = mu * ((Hx + 1/mu*Q )') * U / (I + mu * (U')*U);
end