clc;
clear all;
close all;

Nbi = 120;
x0  = 0;
Nmc = 15000;

a1 = -0.288;
a2 = 0.009;

v     = @(x)a1*x.^2 + a2*x.^4;
p     = @(x)exp(-v(x));
r     = @(x)p(x)*lognrnd(0, sqrt(2));
p_W   = @(w,x)lognpdf(w, 0, sqrt(2));
Xnext = @(x)x+10*(rand(1)-0.5);

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
    p(linspace(h.BinLimits(1),h.BinLimits(2),h.NumBins))/...
    sum(p(linspace(h.BinLimits(1),h.BinLimits(2),h.NumBins))*h.BinWidth),...
    'LineWidth',3,'Color','r');
legend('Pseudo marginal', 'Naive sampling', 'Real distribution');
title('Pseudo marginal vs. Naive sampling');