function [m, NQT, PMF] = heuristic_fec( M, k, tree, PLR, Param )

k = k(:);
N = length(k); % number of frames in the intra-period
m = zeros(N,1);

if PLR == 0
    NQT = ;
    PMF = zeros(1,N+1); PMF(end) = 1;
    return
end

m(1:4:end) = M/
phi = zeros(N,1);
for i = 1:N
    phi(i) = sum(binopdf(0:m(i),m(i)+k(i),PLR));
end

[NQT, PMF] = calculate_mean_nqt( tree, phi(:,1), Param.alpha_t, Param.tmax, Param.ipr );

return