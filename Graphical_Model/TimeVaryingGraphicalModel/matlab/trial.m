% StARS_basis_expansion.m
% Using Stability Selection for choosing tuning paperameter in the
% functional coefficient model through basis expansion with total variation
% penalty.

%% data processing
load mat_data/3dtoyexample1.mat

Data = struct2table(Toydata);
nodename = 'X1';

X = Data(:, 1:(end-1));
X(:, nodename) = [];
X = table2array(X);
Y = table2array(Data(:,nodename));
time = unique(Data.Time);
grid = 0:0.05:1;
set = 0:0.001:1;
%%
W = 0.01:0.02:0.09;
L = [0.1,0.5,1:1:20,21:2:29,30:5:50];
count_result = ones(length(L)*length(W), 3);
support_result = ones(length(L)*length(W), 3);


%%
[opt_bet1, opt_gam_mat1]= basis_expansion_tv( X, Y, time, grid, order, w0, W(1), set, 11 );
[opt_bet2, opt_gam_mat2]= basis_expansion_tv( X, Y, time, grid, order, w0, W(2), set, 10 );
[opt_bet3, opt_gam_mat3]= basis_expansion_tv( X, Y, time, grid, order, w0, W(3), set, 11 );
[opt_bet4, opt_gam_mat4]= basis_expansion_tv( X, Y, time, grid, order, w0, W(4), set, 17 );
%%
Lambda = [0.1,0.5,1:1:20,21:2:29,30:2:50];
count_result = ones(length(Lambda), 3);
support_result = ones(length(Lambda), 3);
for l = 1:length(Lambda)
    [  ~, stab_count_val, ~, stab_support_val] = stars_basis_expansion( Data, nodename, M, Lambda(l), grid, order, 1, 0.05, 0.2 );
    count_result(l,:) = [stab_count_val, Lambda(l),w];
    support_result(l,:) =  [stab_support_val, Lambda(l),w];
end

%%
save mat_data/stars1.mat count_result support_result