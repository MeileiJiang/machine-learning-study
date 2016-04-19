function [] = solve_kernel_group(n, num_rep)
 
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
lambda_max = 3*lambda_max;
while lambda_max > 0.0001
  sol_new.lambda = [sol_new.lambda, lambda_max];
  lambda_max = lambda_max / 1.04;
end
sol_new.lambda = fliplr(sol_new.lambda);
%%%%
sol_new.beta = {};
sol_new.h = h;

it_l = 0;
for lambda=sol_new.lambda
  it_l = it_l + 1;
  fprintf(1, 'lambda=%.4f %d\n', lambda,it_l);
  if it_l==1    
    sol_new.beta{it_l} = sol.beta{I};     
  else
    sol_new.beta{it_l} = sol_new.beta{it_l-1};
  end  
  iter = 0;  stop_criter = 10;
  while (iter<1000 && stop_criter>1e-4)
    iter = iter + 1;
    if mod(iter, 10) == 1,  
      fprintf(1, '%d %.5f\n', iter, stop_criter);
    end
    beta = sol_new.beta{it_l};
    tmp = sqrt(sum( beta.^2, 2 ));
    ind = find( tmp > 1e-8 );
    indc = setdiff(1:p, ind);
    d = diag((lambda ./ avg_norm(ind)) ./ tmp(ind));
    for T=1:n
      t = ((T/n - (1:n)/n)/h).^2;
      t( find(abs(t) >= 1) ) = 1;
      t = sqrt(0.75/h*(1 - t));
      ty = y.*t';
      tx = diag(t) * x(:, ind);    

      sol_new.beta{it_l}(ind, T) = pinv(tx'*tx + d)*tx'*ty;
      sol_new.beta{it_l}(indc, T) = 0;
    end
    stop_criter = max(max(abs(beta - sol_new.beta{it_l})));
  end
end

save(sprintf('output/kernel_group_lasso/kernel_group_lasso_estim_%d_%d', n, num_rep), 'sol_new', '-mat');
quit
