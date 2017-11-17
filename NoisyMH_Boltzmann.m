function [ Sample ] = NoisyMH_Boltzmann( Nbi,x0, r, p_W, Xnext )
%METROPOLISHASTING Summary of this function goes here
%   Detailed explanation goes here
   %    Nmc: total number of samples
   %    p: target probablity distribution
   %    Xnext: the transition function of X' given X
   %    p_W(w,x): the pdf of W|x
p_WX = @(w,x)r(x)*p_W(w,x);   

X = [x0 ; zeros(Nbi,1)];
W = [0 ; zeros(Nbi,1)];
p_WXnew = p_WX(W(1),X(1)); 
for n=2:Nbi+1
   p_WXold = p_WXnew;
   X(n) = Xnext(X(n-1));
   W(n) = lognrnd(0, sqrt(2));
   p_WXnew = p_WX(W(n),X(n));
   A    = min(1, p_WXnew*p_W(W(n-1),X(n-1))/...
       (p_WXold*p_W(W(n),X(n))) ); % Acceptance probability
   if rand(1)>A
       X(n) = X(n-1);
       W(n) = W(n-1);
       p_WXnew = p_WXold;
   end
end
Sample = X(end);
end

