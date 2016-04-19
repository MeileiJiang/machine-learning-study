function [] = solve_bootstrap(n, num_rep)

x = load(sprintf('data/sample_%d/x_%d', n, num_rep));
y = load(sprintf('data/sample_%d/y_%d', n, num_rep));

[n, p] = size(x);

lambda = max(max(abs( 2 * flipud(cumsum( flipud(x).*repmat( flipud(y), 1, p) )))));
rate = 1.03;

Lambda = []; count = 0;
while lambda > 0.001,
  count = count + 1;
  Lambda(count) = lambda;
  lambda = lambda / rate;
end

for num_bootstrap=0:99
  fprintf(1, 'Bootstrap %d\n', num_bootstrap);
  if num_bootstrap == 0
    res = solve_path(y + normrnd(0, 1e-3, n, 1), x, Lambda);
  else
    res = solve_path(y + normrnd(0, 1e-3, n, 1), x, Lambda, res);
  end
  save(sprintf('output/estim_%d_%d_%d', n, num_rep, num_bootstrap), 'res', '-mat');
end

