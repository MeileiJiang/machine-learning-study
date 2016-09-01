function [] = solve_kernel_lasso(n, num_rep)
 
x = load(sprintf('data/sample_%d/x_%d', n, num_rep));
y = load(sprintf('data/sample_%d/y_%d', n, num_rep));

[n, p] = size(x);
bandwidths = (5:5:40) / n;
load(sprintf('output/kernel/kernel_estim_%d_%d', n, num_rep), '-mat');

[tmp, I] = min(sol.objective);
h = bandwidths(I);
avg_norm = sqrt(sum( sol.beta{I}.^2, 2 )/n);

sol_new = struct();

%%%% find lambda
sol_new.lambda = [];
lambda_max = 0;
for i=5:5:100
  T = i;
  t = ((T/n - (1:n)/n)/h).^2;
  t( find(abs(t) >= 1) ) = 1;
  t = sqrt(0.75/h*(1 - t));
  ty = y.*t';
  tx = diag(t) * x;
  tmp = max(abs(2*tx'*ty) .* avg_norm);
  if tmp > lambda_max, lambda_max = tmp; end
end
lambda_max = 2*lambda_max;
while lambda_max > 1e-8
  sol_new.lambda = [sol_new.lambda, lambda_max];
  lambda_max = lambda_max / 1.07;
end
%%%%
sol_new.beta = {};
sol_new.bic = zeros(n, size(sol_new.lambda, 2));
sol_new.h = h;

C = log(n*h)/(n*h);
for T=1:n
  fprintf(1, 'T=%d\n', T);
  t = ((T/n - (1:n)/n)/h).^2;
  t( find(abs(t) >= 1) ) = 1;
  t = sqrt(0.75/h*(1 - t));
  ty = y.*t';
  tx = diag(t) * x;
  
  beta = zeros(p, size(sol_new.lambda, 2)); it = 0;
  b = zeros(p, 1);
  for lambda=sol_new.lambda
    it = it + 1;
    b = shooting_lasso(tx, ty, lambda ./ avg_norm, b);
    beta(:, it) = b;
    % bic = log(RSS) + df * log(nh)/nh
    sol_new.bic(T, it) = nnz(b)*C + log((norm(ty - tx*b, 'fro')^2)/n);
  end
  sol_new.beta{T} = beta;
end

save(sprintf('output/kernel_lasso/kernel_lasso_estim_%d_%d', n, num_rep), 'sol_new', '-mat');
