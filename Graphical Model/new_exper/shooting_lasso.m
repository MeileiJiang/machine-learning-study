function b = shooting_lasso(X, Y, lambda, b)
% Learn l1-regularized l2 loss using shooting: Y = A * X.
% X: covariate, T observations by p dimensions;
% Y: response variable, T observations by d dimensions;
% lambda: regularization parameter;
% b: regression coefficients;

[T, p] = size(X);
[T, d] = size(Y);

XTX = 2 * X' * X;
XTY = 2 * X' * Y;

counter = 0;
b1 = b;
ob = b1 + 10;
while(counter < 1000 && norm(b1 - ob,'fro') > 1e-3)
  counter = counter+1;
  ob = b1;
  if (mod(counter, 2) == 1)
    for j=1:p
      ind = setdiff(1:p, j);
      C = XTX(j, ind) * b1(ind, 1) - XTY(j);
      if (abs(C) < lambda(j))
	b1(j) = 0;	
      elseif (C < 0)
	b1(j,1) = (-C - lambda(j)) / XTX(j, j);
      else
	b1(j,1) = (-C + lambda(j)) / XTX(j, j);
      end	
    end
    aset = find( abs(b1) > eps );
  else
    counter1 = 0;
    ob = ob + 10;
    while (counter1 < 100 && norm(b1 - ob, 'fro') > 1e-3)     
      counter1 = counter1 + 1;
      ob = b1;
      for j=aset'
	ind = setdiff(1:p, j);
	C = XTX(j, ind) * b1(ind, 1) - XTY(j);
	if (abs(C) < lambda(j))
	  b1(j) = 0;	
	elseif (C < 0)
	  b1(j,1) = (-C - lambda(j)) / XTX(j, j);
	else
	  b1(j,1) = (-C + lambda(j)) / XTX(j, j);
	end		
      end
    end
    ob = b1 + 10;
  end
end

b = b1;