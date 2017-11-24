true_mean = 0;
true_sigma = 1;
% likelihood_func = @(x, mean, sigma) normpdf(x, mean, sigma);
% the above function to calcalate in matrix form, for speed
likelihood_func = @(x, mean, sigma)prod(normpdf(repmat(x,[1 numel(mean)]),repmat(mean, [1 numel(x)])',repmat(sigma,[1 numel(x)])' ), 1);

%% generate data
N=20;
observed_data = normrnd(true_mean, true_sigma, [N 1]);

%% Do importance sampling
% create many samples for mean and sigma
N = 10^6;
mean_samples = (rand(N,1)-0.5)*5; 
sigma_samples = rand(N, 1) * 10; 
proposal = 1/N;

% evaluate likelihood for each (mean, sigma) sample 
target = likelihood_func(observed_data, mean_samples, sigma_samples);

% calculate importance weight
w = target ./ proposal;
w = w ./ sum(w);

% resample, with replacement, according to importance weight
sample_ind = randsample([1:N],N,true,w);

mean_samples = mean_samples(sample_ind);
sigma_samples = sigma_samples(sample_ind);

%% plot
plot(mean_samples,sigma_samples,'.')
xlabel('mean')
ylabel('sigma')
hold on
plot(true_mean, true_sigma, 'r.','MarkerSize',5^2)
axis square