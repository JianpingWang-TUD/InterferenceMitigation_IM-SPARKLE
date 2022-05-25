function x = update_LowRank(y,i,p,beta_1,mu,U,V,Q)
    F1 = U*(V') - 1/mu*Q;
    j = 0;
    F1_flip = fliplr(F1);
    f1 = zeros(size(F1,1)+size(F1,2)-1,1);
   for m = size(F1,2)-1:-1:-(size(F1,1)-1)
       j = j+1;
       f1(j) = mean(diag(F1_flip,m));
   end
%     f1 = [F1(:,1);F1(end,2:end).'];
    x  = 1/(beta_1+mu) * ( beta_1*(y-i+1/beta_1*p) + mu*f1 ); 
end