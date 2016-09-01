function [new_beta] = solve(y, x, lambda_1, lambda_tv, old_solution)

[n, p] = size(x);

converge_eps = 10;
new_beta = old_solution;

iter = 0;

while (iter < 50 && converge_eps > 1e-3)
    iter = iter+1;
    disp(['Iteration ', num2str(iter)]);
    old_beta = new_beta;
    
    for i=1:p
        ind = setdiff(1:p, i);
        residuals = y - sum(x(:, ind) .* new_beta(:, ind), 2);
        cvx_begin
            cvx_precision low
            cvx_quiet(true);
            variable var_beta(n, 1)
       
            if lambda_1 ~= 0
                minimize sum(power(residuals - x(:,i).*var_beta, 2)) + lambda_tv * norm( var_beta(2:end, 1) - var_beta(1:end-1, 1), 1 ) + lambda_1 * norm( var_beta, 1 )
            else
                minimize sum(power(residuals - x(:,i).*var_beta,2)) + lambda_tv ...
                    * norm( var_beta(2:end, 1) - var_beta(1:end-1, 1), 1 ...
                            ) + lambda_tv * abs(var_beta(1,1))
            end
        cvx_end 
        new_beta(:, i) = var_beta;
        new_beta( find( abs(new_beta) < 1e-4 ) ) = 0; 
    end
    converge_eps = max(max(abs( (new_beta - old_beta) )));
    fprintf(1, 'convergence  ==> %.4f\n', converge_eps);
end

