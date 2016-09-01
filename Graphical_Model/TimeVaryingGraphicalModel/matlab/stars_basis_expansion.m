function [ count, stab_count_val, support, stab_support_val] = stars_basis_expansion( Data, nodename, M, lambda, grid, order, w0, w, frac)
% Using Stability Selection for choosing tuning paperameter in the
% functional coefficient model through basis expansion with total variation
% penalty.

time = unique(Data.Time);
set = 0:0.001:1;

count = [];
support = [];
for k = 1:M
    subdata = subsample(Data, frac);
    X = subdata(:, 1:(end-1));
    X(:,nodename) = [];
    X = table2array(X);
    Y = table2array(subdata(:,nodename));
    
    [bet, gam_mat] = basis_expansion_tv( X, Y, time, grid, order, w0, w, set, lambda );
    if isempty(count)
        count = abs(gam_mat) > 0;
    else
        count = count +  (abs(gam_mat) > 0);
    end
    
    if isempty(support)
        support = abs(bet) > 0;
    else
        support = support +  (abs(bet) > 0);
    end
end

count = count/M;
support = support/M;

stab_count = 2*count.*(1-count);
stab_count_val = mean2(stab_count);
 
stab_support = 2*support.*(1-support);
stab_support_val = mean2(stab_support);
end

