function [b] = solve_tv_shooting(y, x, lambda, b0)
% this version is slover than solve_tv.m


[n, p] = size(x);

if (nargin == 3),
  XTX = 2 * x' * x;
  XTY = 2 * x' * y;
  b = (XTX + 2*lambda*eye(p)) \ XTY;
  b = repmat(b, 1, n);
else
  b = b0;
end

XTY = 2 * flipud(cumsum( flipud(x).*repmat( flipud(y), 1, p) ));
XTX = zeros(n, p);

counter = -1;
ob = b + 10;
tmp = 0; tmp1 = 0;
while(counter < 1000 && norm(b - ob,'fro') > 1e-4)
  counter = counter + 1;
  fprintf(1, 'Counter: %d Err %.5f', counter, norm(b - ob,'fro'));
  ob = b;
  if (mod(counter, 2) == 0)
    % do full update
    fprintf(1, ' full\n');
    for k=1:p
      kset = setdiff(1:p, k);
      XTX(:, kset) = x(:, kset) .* b(kset, :)' .* repmat(x(:, k), 1, p - 1);
      XTX(:, k) = zeros(n, 1);
      XTX = 2 * flipud( cumsum( flipud(XTX) ) );
      sXTX = sum( XTX, 2 );
      b_prime = b(k, :) - [0, b(k, 1:end-1)];     
      for l=1:n     
	tmp = 0;
	for l_p=l:n
	  lset = setdiff(1:l_p, l);
	  tmp = tmp + x(l_p, k)^2 * sum(b_prime(lset));
	end
	tmp = 2 * tmp;
	tmp = tmp + sXTX(l) - XTY(l, k);
	
	if (abs(tmp) <= lambda)
	  b_prime(l) = 0;
	elseif (tmp < 0)
	  b_prime(l) = ( -tmp - lambda) / (2* sum(x(l:end, k).^2));
	else
	  b_prime(l) = ( -tmp + lambda) / (2* sum(x(l:end, k).^2));
	end
      end
      b(k, :) = cumsum(b_prime);
    end
    [aset_r, aset_c] = find(abs(b) > eps);
  else
    % do active set update
    b1 = b + 10;
    counter1 = 0;
    num_elem = size(aset_r, 1);
    fprintf(1, ' active %d\n',num_elem);
    while(counter1 < 100 && norm(b - b1,'fro') > 1e-4)    
      counter1 = counter1 + 1;
      k = -1;
      b1 = b;
      for i=1:num_elem
	if (k ~= aset_r(i))
	  k = aset_r(i);
	  kset = setdiff(1:p, k);
	  XTX(:, kset) = x(:, kset) .* b(kset, :)' .* repmat(x(:, k), 1, p - 1);
	  XTX(:, k) = zeros(n, 1);
	  XTX = 2 * flipud( cumsum( flipud(XTX) ) );
	  sXTX = sum( XTX, 2 );
	  b_prime = b(k, :) - [0, b(k, 1:end-1)];     	
	end
	l = aset_c(i);
	tmp = 0;
	for l_p=l:n
	  lset = setdiff(1:l_p, l);
	  tmp = tmp + x(l_p, k)^2 * sum(b_prime(lset));
	end
	tmp = 2 * tmp;
	tmp = tmp + sXTX(l) - XTY(l, k);
	
	if (abs(tmp) <= lambda)
	  b_prime(l) = 0;
	elseif (tmp < 0)
	  b_prime(l) = (-tmp - lambda) / (2*sum(x(l:end, k).^2));
	else
	  b_prime(l) = (-tmp + lambda) / (2*sum(x(l:end, k).^2));
	end      
	b(k, :) = cumsum(b_prime);
      end
    end
  end
end
