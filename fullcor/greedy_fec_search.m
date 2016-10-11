function [best_alloc, best_nqt] = greedy_fec_search( M, k, Markov, Param, varargin )

k = k(:);
N = length(k); % number of frames in the intra-period

if ~isempty(varargin)
    in = varargin{1}; in = in(:);
else
    in = zeros(length(k),1);
end

best_alloc = in;
% best_nqt = monte_carlo( markov, k, best_alloc, tree, param, 2000 );
% best_nqt = mnqt( (0:N)/param.ipr, param.alpha_t, param.tmax ) * calculate_pmf( tree, k, best_alloc, markov, markov.ss, 1, N, [] );
best_nqt = compute_mean_nqt( Param, Markov, k, best_alloc );
while M > 0
    
    curr_alloc = best_alloc;
    seen = cell(1,Param.L);
    for j = 1:N
        if any(seen{Param.layermap(j)} == k(j))
            continue;
        else
            seen{Param.layermap(j)} = [seen{Param.layermap(j)},k(j)]; 
        end
        m = curr_alloc;
        m(j) = m(j) + 1;

%         nqt = mnqt( (0:N)/param.ipr, param.alpha_t, param.tmax ) * calculate_pmf( tree, k, m, markov, markov.ss, 1, N, [] );
%         nqt = monte_carlo( markov, k, m, tree, param, 2000 );
        nqt = compute_mean_nqt( Param, Markov, k, m );
        if nqt >= best_nqt
            best_nqt = nqt;
            picked = j;
        end        
    end
    
    best_alloc(picked) = best_alloc(picked)+1;
    
    M = M-1;
end

return