function [x, intf,  rerr] = lowRaS_Hankel(y, R_mtr, beta_1, mu,...
                            tau, k_beta, k_mu)
% lowRaS_Hankel--- implements the signal separation by decomposing a Hankel  
%                   matrix as the sum of one low-rank and one sparse Hankel 
%                   matrices
%
% Input Parameters:
%     y --- the signal vector, column vector
%     R_mtr --- the approximate of the Hankel matrix formed by usable signal
%     beta_1 --- 
%     mu ---
%     tau ---
%
% Output:
%     x --- the usable signal obtained after separation
%    intf --- the interference
%
% $$
% Created by Jianping Wang @ MS3, TU Delft, May 15, 2020
% email: J.Wang-4@tudelft.nl or jianpingwang87@gmail.com
%==========================================================================
%

y_len = length(y);
if size(y,1)==1
    y = y.';
end


x = zeros(y_len,1);
intf = zeros(y_len,1);

p = zeros(y_len,1);
Q = zeros(floor(y_len/2), y_len-floor(y_len/2)+1);
% R_mtr = min(floor(length(x))/2, length(x)-floor(length(x)/2)+1);
               
U = randn(floor(y_len/2), R_mtr);
V = randn(y_len-floor(y_len/2)+1, R_mtr);
% Hy = hankel(y(1:floor(length(y)/2)),y(floor(length(y)/2):end));
% m = 0;
% while norm(U*V'-Hy,'fro')/norm(Hy,'fro') >1e-7
%     m = m+1;
%     fprintf('%d %12.8d\n', m, norm(U*V'-Hy,'fro')/norm(Hy,'fro'))
%     U = Hy*V/(V'*V);
%     V = Hy'*U/(U'*U);
% end
Hx = hankel(x(1:floor(y_len/2)),x(floor(y_len/2):end));
%% ======== Update =========
rerr = zeros(1,1000);

k = 0;
while norm(y-x-intf)/norm(y)>1e-6 && k<1000%1e-8
    k = k+1;
    if mod(k,10)==0
        fprintf('%d %d %15.8d\n', k, rank(Hx), norm(y-x-intf)/norm(y))
        beta_1 = beta_1*k_beta; %1.6   1.3
    end
    rerr(k) = norm(y-x-intf)/norm(y);
    
   % Update x
    x = update_LowRank(y,intf,p,beta_1,mu,U,V,Q);
   % Update i --- interference
    intf = update_In(y,x,tau,p,beta_1);
   % Update U
    [U,Hx] = update_U(x,V,Q,mu);
   % Update V
    V = update_V(Hx,U,Q,mu);
   % Update p
    p = p + beta_1*(y-x-intf);
   % Update Q
    Q = Q + mu*(Hx-U*(V'));
%     beta_1 = beta_1*1.05;
    mu = mu*k_mu;%1.2; %1.05
end

rerr = rerr(1:k);