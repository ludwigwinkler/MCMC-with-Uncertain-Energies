
nb_samples = 1000000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

variance = 2;
N = 128;
samples = metropolis(nb_samples, variance, N);
plot_boltzmann(1, samples(:,1), 'Metropolis-Hastings: Boltzmann Distribution')
plot_boltzmann(2, samples(:,2), 'Metropolis-Hastings (case 1): $\sigma^2 = 2$')
plot_boltzmann(3, samples(:,3), 'Metropolis-Hastings (case 2): $\sigma^2 = 2$')
plot_boltzmann(4, samples(:,4), 'Naive Metropolis-Hastings: $\sigma^2 = 2$')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_boltzmann(fig_nb, samples, title_txt)
figure(fig_nb)
start = 2000;
edges = linspace(min(samples),max(samples),100);
counts = histcounts(samples(start:end), edges);
histogram('BinEdges',edges, 'BinCounts', counts/sum(counts));
xlabel('state $x$', 'Interpreter','LaTex')
ylabel('target density $\pi(x)$', 'Interpreter','LaTex')
title(title_txt, 'Interpreter','LaTex')
hold on
x = linspace(min(samples),max(samples),100);
y = exp(-arrayfun(@V,x));
plot(x, y/sum(y), 'Color', 'red', 'LineWidth',2);
hold off
end

function x = metropolis(nb_samples, variance, N)
x = zeros(nb_samples,4);
for i=2:nb_samples
    
    u = random('Uniform', -0.5, 0.5);
    u_accept = random('Uniform', 0, 1);
    
    %%% without noise %%%
    x_prev = x(i-1,1);
    x_star = x_prev + u; 
    Delta = V(x_star) - V(x_prev);
    q = exp(-Delta);
    x(i,1) = helper(q, x_prev, x_star, u_accept);
    
    %%% with noise naive %%%
    x_prev = x(i-1,4);
    x_star = x_prev + u; 
    Delta = V(x_star) - V(x_prev) + random('Normal', 0, sqrt(variance));
    q = exp(-Delta);
    x(i,4) = helper(q, x_prev, x_star, u_accept);
    
    %%% with noise case 1 %%%
    x_prev = x(i-1,2);
    x_star = x_prev + u; 
    Delta = V(x_star) - V(x_prev);
    delta = random('Normal', Delta, sqrt(variance));
    q = exp(-delta-variance/2);
    x(i,2) = helper(q, x_prev, x_star, u_accept);
    
    %%% with noise case 2 %%%
    x_prev = x(i-1,3);
    x_star = x_prev + u; 
    Delta = V(x_star) - V(x_prev);
    samples_delta = random('Normal', Delta, sqrt(N * variance), [N,1]);
    
    
    Delta_est = mean(samples_delta);
    variance_est = 1/(N*(N-1)) * sum( (samples_delta - Delta_est).^2 );
    
    %q = exp(-Delta_est-variance_est/2);

    b = bessel(N, variance_est);
    q = exp(-Delta_est-b);
    x(i,3) = helper(q, x_prev, x_star, u_accept);
end
end


function x = helper(q, x_prev, x_star, u_accept)

if q > 1
    A = 1;
else
    A = q;
end

if u_accept < A
    x = x_star;
else
    x = x_prev;
end

end


function b = bessel(N, variance)
b = variance/2 + variance^2/(4*(N+1)) + variance^3/( 3*(N+1)*(N+3)  );
end

function v = V(x)
a1 = -0.288;
a2 = 0.009;
v = a1 * x^2 + a2 * x^4;
end