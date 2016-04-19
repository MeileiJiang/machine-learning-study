function [lambda, cum_freq_br] = plot_break(n, num_rep)

lambda = [];
threshold = 2;

for n_boot=0:99
  load(sprintf('output/estim_%d_%d_%d', n, num_rep, n_boot), '-mat');
  if n_boot==0
    lambda = res.lambda;
    n_lambda = size(lambda, 2);
    cum_freq_br = zeros(n, n_lambda);
  end
  for i=1:n_lambda
    [freq_br, delta] = calculate_breaks(res.beta{i}, 5e-3);
    ind = find( freq_br > threshold);
    freq_br( ind ) = 1;
    freq_br( setdiff(1:n, ind) ) = 0;
    cum_freq_br(:, i) = cum_freq_br(:, i) + freq_br';
  end  
end

plot(1:n_lambda, cum_freq_br);