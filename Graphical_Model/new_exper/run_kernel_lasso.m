function [] = run_kernel_lasso(n)

  for rep=0:19
    fprintf(1, 'rep=%d\n', n);
    solve_kernel_lasso(n, rep);
  end
