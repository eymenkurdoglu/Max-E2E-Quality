function [best_alloc, best_nqt, best_pmf] = greedy_fec_search( M, k, tree, PLR, Param, varargin )

k = k(:);
N = length(k); % number of frames in the intra-period

if ~isempty(varargin)
    in = varargin{1}; in = in(:);
else
    in = zeros(length(k),1);
end

if PLR == 0
    best_alloc = in;
    best_nqt = N;
    best_pmf = zeros(1,N+1); best_pmf(end) = 1;
    return
end

phi = zeros(N,2);
for i = 1:N
    % First column : frame arrival probabilities
    phi(i,1) = sum(binopdf( 0:in(i)  , in(i)+k(i)  , PLR));
    % Second column: frame arrival probabilities with an additional FEC block
    phi(i,2) = sum(binopdf( 0:in(i)+1, in(i)+k(i)+1, PLR));
end

best_alloc = in;
[best_nqt, best_pmf] = calculate_mean_nqt( tree, phi(:,1), Param.alpha_t, Param.tmax, Param.ipr );

while M > 0
    % check each neighbor
    for j = randperm(N)
        
        phi_    = phi(:,1);
        phi_(j) = phi(j,2);
        
        [nqt, pmf] = calculate_mean_nqt( tree, phi_, Param.alpha_t, Param.tmax, Param.ipr );
        
        if nqt >= best_nqt
            best_pmf = pmf;
            best_nqt = nqt;
            picked = j;
        end
    end
    
    best_alloc(picked) = best_alloc(picked)+1;
    
    % need to only change the success prob corr. to the picked frame
    phi(picked,1) = phi(picked,2);
    
    % to update the value in the second column
    phi(picked,2) = sum( binopdf(0:(best_alloc(picked)+1), k(picked)+best_alloc(picked)+1, PLR) );
    
    M = M-1;
end

return