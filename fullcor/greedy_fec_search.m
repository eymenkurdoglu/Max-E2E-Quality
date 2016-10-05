function [best_alloc, best_nqt] = greedy_fec_search( f, k, tree, eps, alpha_t, tmax, ipr, varargin )
% GREEDY_FEC_SEARCH    greedy search for distributing f FEC packets on
% frames of length k (should be faster than GREEDY_FEC_SEARCH3 for independent losses)
% f: total number of FEC blocks we can use
% k: vector of frame lengths in the intra-period 
% tree: frame dependency tree struct
% eps: iid packet loss rate
% alpha_t, tmax: parameters of NQT
% ipr: intra-period in seconds
% in: initial FEC block distribution (varargin)
% =====================================================================

k = k(:);
N = length(k); % number of frames in the intra-period

if ~isempty(varargin)
    in = varargin{1}; in = in(:);
else
    in = zeros(length(k),1);
end

if eps == 0
    best_alloc = in;
    best_nqt = N;
    return
end

phi = zeros(N,2);
for i = 1:N
    % First column : frame arrival probabilities
    phi(i,1) = sum(binopdf( 0:in(i)  , in(i)+k(i)  , eps));
    % Second column: frame arrival probabilities with an additional FEC block
    phi(i,2) = sum(binopdf( 0:in(i)+1, in(i)+k(i)+1, eps));
end

best_alloc = in;
best_nqt = calculate_mean_nqt( tree, phi(:,1), alpha_t, tmax, ipr );

while f > 0
    % check each neighbor
    for j = randperm(N)
        
        phi_    = phi(:,1);
        phi_(j) = phi(j,2);
        
        nqt = calculate_mean_nqt( tree, phi_, alpha_t, tmax, ipr );
        
        if nqt >= best_nqt
            best_nqt = nqt;
            picked = j;
        end
    end
    
    best_alloc(picked) = best_alloc(picked)+1;
    
    % need to only change the success prob corr. to the picked frame
    phi(picked,1) = phi(picked,2);
    
    % to update the value in the second column
    phi(picked,2) = sum( binopdf(0:(best_alloc(picked)+1), k(picked)+best_alloc(picked)+1, eps) );
    
    f = f-1;
end

return