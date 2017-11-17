clc;
clear all;
close all;

Nbi = 120;
x0  = 0;
Nmc = 15000;

p     = @(x)1/sqrt(2*pi)*exp(-x.^2/2);
r     = @(x)p(x)*(1+0.9*(exprnd(1)-1));
p_W   = @(w,x)1/0.9*exp(-(w-0.1)/0.9)*double(w>=0.1);
Xnext = @(x)x+rand(1)-0.5;

Samples  = zeros(Nmc,1);
NoisySam = zeros(Nmc,1);
for n=1:Nmc
    NoisySam(n) = NoisyMH( Nbi,x0, r, p_W, Xnext );
    Samples(n)  = metropolishasting( Nbi,x0, r, Xnext );
end

Nh = histogram(NoisySam,'NumBins',50,'Normalization','pdf');
hold on;
h = histogram(Samples,'NumBins',50,'Normalization','pdf');
plot(linspace(h.BinLimits(1),h.BinLimits(2),h.NumBins),...
    p(linspace(h.BinLimits(1),h.BinLimits(2),h.NumBins)),'LineWidth',3,'Color','r');
legend('Pseudo marginal', 'Naive sampling', 'Real distribution');
title('Pseudo marginal vs. Naive sampling');



