function i = update_In(y,x,tau,p,beta_1)
    f = y-x+1/beta_1*p;
    i = soft(f,tau/beta_1);
end