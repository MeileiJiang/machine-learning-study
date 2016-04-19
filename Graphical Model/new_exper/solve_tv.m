function [beta] = solve_tv(y, x, lambda, init_x)

eps = 1e-3;
  
[n, p] = size(x);
%p = p + 1;
%x = [ones(n, 1) / sqrt(n) x];

beta_tilde = init_x - [zeros([1, p]); init_x(1:end-1, :)];

% precalculate
x_ssum = zeros([n, p]);
for i=1:n   
  x_ssum(i, :) = sum( x(i:end, :) .* x(i:end, :), 1 );
end

% start
beta = cumsum(beta_tilde, 1);
iteration = 0;
while (iteration < 50)
  iteration = iteration + 1;
  new_beta_tilde = beta_tilde;  

  if (mod(iteration, 2) == 1)
    for j=1:p
      disp([sprintf('Iteration %d FULL feature %d', iteration, j)]);
      ind = (1:p);
      ind = ind(1:p ~= j);        
      y_tilde = y - sum(x(:, ind) .* beta(:, ind), 2);     
      
      new_beta_tilde(:, j) = compute_column_m(y_tilde, x(:, j), new_beta_tilde(:, j), x_ssum(:, j), lambda, eps);
      beta(:, j) = cumsum(new_beta_tilde(:, j));
    end
    active_set = find( sum( abs(beta), 1 ) ~= 0 );
    if max(max( abs(new_beta_tilde - beta_tilde) )) > eps
      beta_tilde = new_beta_tilde;
    else
      break;
    end  
  else
    count = 0;
    while (count < 100)
      count = count + 1;
      for j=active_set
	disp([sprintf('Iteration %d %d ACTIVE feature %d', iteration, count, j)]);
	ind = (1:p);
	ind = ind(1:p ~= j);        
	y_tilde = y - sum(x(:, ind) .* beta(:, ind), 2);     
      
	new_beta_tilde(:, j) = compute_column_m(y_tilde, x(:, j), new_beta_tilde(:, j), x_ssum(:, j), lambda, eps);
	beta(:, j) = cumsum(new_beta_tilde(:, j));
      end
      if max(max( abs(new_beta_tilde - beta_tilde) )) > eps
	beta_tilde = new_beta_tilde;
      else
	break;
      end  
    end
  end  
end

beta = beta;
