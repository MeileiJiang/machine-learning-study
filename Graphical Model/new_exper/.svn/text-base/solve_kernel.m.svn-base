function [] = solve_kernel(n, num_rep)

x = load(sprintf('data/sample_%d/x_%d', n, num_rep));
y = load(sprintf('data/sample_%d/y_%d', n, num_rep));

[n, p] = size(x);
bandwidths = (5:5:40) / n;
num_band = size(bandwidths, 2);
sol = struct();
sol.n = n;
sol.p = p;
sol.rss = zeros(num_band, n);
sol.df = zeros(num_band, n);
sol.beta = {};
sol.objective = zeros(1, num_band);

it_h = 0;
for h=bandwidths
  it_h = it_h + 1;
  beta = zeros(p, n);
  for T=1:n
    t = ((T/n - (1:n)/n)/h).^2;
    t( find(abs(t) >= 1) ) = 1;
    t = sqrt(0.75/h*(1 - t));
    ty = y.*t';
    tx = diag(t) * x;
    tmp = pinv(tx'*tx)*tx';
    beta(:, T) = tmp*y;
    sol.rss(it_h, T) = norm(y - tx*beta(:, T), 'fro')^2;
    sol.df(it_h, T) = trace(tx*tmp);
  end
  sol.beta{it_h} = beta;
  sol.objective(1, it_h) = sum(sol.rss(it_h, :) ./ (n*(1-sol.df(it_h, :)/n).^2));
end

save(sprintf('output/kernel/kernel_estim_%d_%d', n, num_rep), 'sol', '-mat');