function [ Sample ] = metropolishasting( Nbi,x0, p, Xnext )
%METROPOLISHASTING Summary of this function goes here
%   Detailed explanation goes here
   %    Nmc: total number of samples
   %    p: target probablity distribution
   %    Xnext: the transition function of X' given X
   
X = [x0 ; zeros(Nbi,1)];
for n=2:Nbi+1
   X(n) = Xnext(X(n-1));
   A    = min(1, p(X(n))/p(X(n-1)) ); % Acceptance probability
   if rand(1)>A
       X(n) = X(n-1);
   end
end
Sample = X(end);
end

