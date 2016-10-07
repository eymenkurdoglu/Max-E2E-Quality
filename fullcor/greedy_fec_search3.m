function [best_alloc, best_nqt] = greedy_fec_search3( f, k, tree, markov, param, varargin )

k = k(:);
N = length(k); % number of frames in the intra-period

if ~isempty(varargin)
    in = varargin{1}; in = in(:);
else
    in = zeros(length(k),1);
end

best_alloc = in;
tic
best_nqt = monte_carlo( markov, k, best_alloc, tree, param, 2000 );
% best_nqt = mnqt( (0:N)/param.ipr, param.alpha_t, param.tmax ) * calculate_pmf( tree, k, best_alloc, markov, markov.ss, 1, N, [] );
toc
while f > 0
    
    curr_alloc = best_alloc;
    
    for j = randperm(N)
        
        m = curr_alloc;
        m(j) = m(j) + 1;

%         nqt = mnqt( (0:N)/param.ipr, param.alpha_t, param.tmax ) * calculate_pmf( tree, k, m, markov, markov.ss, 1, N, [] );
        tic
        nqt = monte_carlo( markov, k, m, tree, param, 2000 );
        toc
        if nqt >= best_nqt
            best_nqt = nqt;
            picked = j;
        end        
    end
    
    best_alloc(picked) = best_alloc(picked)+1;
    
    f = f-1;
end

return