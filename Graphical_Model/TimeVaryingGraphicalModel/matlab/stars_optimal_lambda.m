function [ opt_lambda ] = stars_optimal_lambda( Data, nodename, M, order, w, grid, Lambda, frac )
%Output opt_lambda = inf_{lambda \in Lambda} sup_{lam > lambda} stars(lam)
%< 0.05 for the function coefficient model with order th degree B-spline
%basis over the knots set grid. The propertion of total variation penalty
%over l1 norm penalty is w.

count_result = ones(length(Lambda), 3);
support_result = ones(length(Lambda), 3);

tic
parfor l = 1:length(Lambda)
    [  ~, stab_count_val, ~, stab_support_val] = stars_basis_expansion( Data, nodename, M, Lambda(l), grid, order, 1, w, frac );
    count_result(l,:) = [stab_count_val, Lambda(l),w];
    support_result(l,:) =  [stab_support_val, Lambda(l),w];
end
toc


end

