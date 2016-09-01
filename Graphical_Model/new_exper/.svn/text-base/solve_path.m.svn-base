function [res] = solve_path(y, x, Lambda, res_old)

res = struct();
res.lambda = Lambda;
res.beta = {};

[n,p] = size(x);

beta = zeros(n, p);
count = 0;
tTotal = tic;
for lambda=Lambda
  count = count + 1;
  tInner = tic;
  if nargin==3
    beta = solve_tv(y, x, lambda, beta);
  else
    beta = solve_tv(y, x, lambda, res_old.beta{count});
  end
  tElapsedInner = toc(tInner);
  fprintf(1, 'Time elapsed: %.3f tElapsed lambda: %.3f\n', tElapsedInner, lambda);
  res.beta{count} = beta;
end
tElapsed = toc(tTotal);
fprintf(1, 'Total time: %.3f\n', tElapsed);