function [new_beta_tilde] = compute_column_m(y_tilde, x, beta_tilde, x_ssum, lambda, eps)

n = length(y_tilde);

new_beta_tilde = beta_tilde;

while true
  for k=1:n
    tmp = cumsum(new_beta_tilde);
    tmp = tmp(k:n, 1) - new_beta_tilde(k);
    y_ttilde = y_tilde(k:n, 1) - x(k:n, 1).*tmp;
    new_beta_tilde(k) = soft_thr( (y_ttilde'*x(k:n, 1)) / x_ssum(k, 1), lambda / x_ssum(k, 1) );
  end
  if max( abs(new_beta_tilde - beta_tilde) ) > eps
    beta_tilde = new_beta_tilde;
  else
    break;
  end
end
