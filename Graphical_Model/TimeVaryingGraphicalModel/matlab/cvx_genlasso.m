% Generalized lasso by cvx
load /Users/meileijiang/researchspace/Graphical_Model/TimeVaryingGraphicalModel/R/3d_Case/temp_rdata/Genlasso_cvx.mat
%%
n = 38;
lambda = 7;
% cvx
cvx_begin
    variable x(n)
    minimize( 0.5*sum_square(U*x - Y) + lambda*norm(A*x,1))
cvx_end
x
plot(A*x)
%%
save /Users/meileijiang/researchspace/Graphical_Model/TimeVaryingGraphicalModel/R/3d_Case/temp_rdata/cvx_result.mat x