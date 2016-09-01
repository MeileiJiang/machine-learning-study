function [] = gen_data(n,p,njumps,dir)

sigma = 0.1;
out_dir = sprintf('data/%s/', dir);
mkdir(out_dir);

params = struct();
params.n = n;
params.p = p;
params.njumps = njumps;
params.s = 5;
params.beta = zeros(p, n);


sign_m = [-1; 1];
block_len = n / njumps;
s = 5;

for bl=1:njumps
  ind = randperm(p);
  ind = ind(1:s);
  val = unifrnd(0.1, 1, s, 1);
  val = val .* randsample(sign_m, s, true);
  params.beta(ind, (block_len*(bl-1)+1+block_len/2) : (block_len*bl)) = ...
      repmat(val, 1, block_len / 2);     
end

save( sprintf('%s/model', out_dir), 'params');

sigma_n = sigma / sqrt(n);
for i=0:49
  x = normrnd(0, 1, n, p);
  err = normrnd(0, sigma_n, n, 1);
  y = sum(x.*params.beta', 2) + err;
  
  save(sprintf('%s/x_%d', out_dir, i), 'x', '-ASCII');
  save(sprintf('%s/y_%d', out_dir, i), 'y', '-ASCII');
end