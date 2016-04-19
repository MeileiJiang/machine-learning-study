function [] = run_tv(n, num_rep)

x = load(sprintf('data/sample_%d/x_%d', n, num_rep));
y = load(sprintf('data/sample_%d/y_%d', n, num_rep));

[n, p] = size(x);

beta = zeros(n, p);
grid_l1 = [];
grid_ltv = [];

l1 = 0.4; it_l1 = 0; 
while l1 > 0.05
  it_l1 = it_l1 + 1;
  ltv = 2.0; it_ltv = 0;
  grid_l1(it_l1) = l1;
  while ltv > 0.6
    it_ltv = it_ltv + 1;
    grid_ltv(it_ltv) = ltv;
    beta = solve(y, x, l1, ltv, beta);
    save(sprintf('output/tv/estim_tv_%d_%d_%d_%d', n, num_rep, it_l1, it_ltv), 'beta', 'l1', 'ltv', '-mat');
    ltv = ltv / 1.2;
  end
  l1 = l1 / 1.15;
end

save(sprintf('output/tv/grid_%d_%d', n, num_rep), 'grid_l1', 'grid_ltv', '-mat');
quit;
