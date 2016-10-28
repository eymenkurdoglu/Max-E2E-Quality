function [best_alloc, best_nqt, best_pmf] = greedy_fec_search( M, k, Markov, Param, varargin )

k = k(:);
N = length(k); % number of frames in the intra-period

if ~isempty(varargin)
    in = varargin{1}; in = in(:);
else
    in = zeros(length(k),1);
end

best_alloc = in;
% best_nqt = monte_carlo( markov, k, best_alloc, tree, param, 2000 );
[best_nqt, best_pmf] = compute_mean_nqt( Param, Markov, k, best_alloc );
while M > 0
    curr_alloc = best_alloc;
    checked = cell(1,Param.L);
    for j = 1:N
        m = curr_alloc;
        m(j) = m(j) + 1;
        
        if ~isempty( checked{Param.layermap(j)} ) && ...
                any( checked{Param.layermap(j)}(1,:) == k(j) & checked{Param.layermap(j)}(2,:) == m(j) )
            continue;
        else
            checked{Param.layermap(j)} = [checked{Param.layermap(j)},[k(j);m(j)]]; 
        end

%         nqt = monte_carlo( markov, k, m, tree, param, 2000 );
        [nqt, pmf] = compute_mean_nqt( Param, Markov, k, m );
        if nqt >= best_nqt
            best_pmf = pmf;
            best_nqt = nqt;
            picked = j;
        end        
    end
    
    best_alloc(picked) = best_alloc(picked)+1;
    
    M = M-1;
end

return