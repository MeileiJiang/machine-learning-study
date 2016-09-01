function [ bet, gam_mat ] = basis_expansion_tv( X, Y, time, grid, order, w0, w, set, lambda )
%basis_expansion_tv solves functional coefficent model through basis
%expansion approach and adds total variation panelty on coefficent
%functions.
%   Coefficient function is approximated by B-spline basis.

% bspline basis matrix
B = bspline_basismatrix(order, grid, time);

% basis expansion of varying coefficient model
n = size(X, 1)/length(time);
p = size(X, 2);

U = zeros(size(X, 1), size(X, 2)*size(B, 2));
for t = 1:length(time)
    U((1+(t-1)*n) : (t*n),:) = X((1+(t-1)*n) : (t*n),:) * kron(eye(p), B(t,:));
end

% penalty matrix
if order > 1
    delta = grid(2) - grid(1);
    b1 = bspline_basismatrix(1, grid, time);
    db = diff(b1.', order-1).'/delta^(order-1);
    DB = diff(db);
else
    DB = diff(B);
end
A0 = [w0*B;w*DB];
A = kron(eye(p), A0);

% cvx to solve genralized lasso problem
cvx_begin quiet
    variable gam(size(A,2))
    minimize( 0.5*sum_square(U*gam - Y) + lambda*norm(A*gam,1))
cvx_end

% threshold the extreme small values in gam
th = 10^(-4);
gam(abs(gam) < th) = 0;
gam_mat = reshape(gam, [],p);

% value of beta
basis = bspline_basismatrix(order, grid, set);
bet = basis*gam_mat;

end

