function [freq_break, delta] = calculate_breaks(beta, eps)

if nargin==1
  eps = 1e-3;
end

[n, p] = size(beta);
delta = beta - [zeros(1, p); beta(1:end-1, :)];
freq_break = delta;

ind = find( abs(delta) > eps );
indc = find( abs(delta) <= eps );

freq_break( ind ) = 1;
freq_break( indc ) = 0;

freq_break = sum(freq_break, 2)';