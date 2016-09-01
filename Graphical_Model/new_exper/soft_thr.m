function b = soft_thr(beta_hat, lambda)
  
if lambda >= abs(beta_hat)
  b = 0;
else
  if beta_hat > 0
    b = beta_hat - lambda;
  else
    b = beta_hat + lambda;
  end;
end