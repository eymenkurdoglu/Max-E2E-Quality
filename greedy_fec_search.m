function [best_alloc, best_score] = greedy_fec_search( f, k, in, templayers, eps )
% GREEDY_FEC_SEARCH    greedy search for distributing f FEC packets on
% frames of length k (should be faster than GREEDY_FEC_SEARCH3 for independent losses)
% f: total number of FEC blocks we can use
% k: vector of frame lengths in the intra-period 
% in: initial FEC block distribution
% templayers: number of temporal layers of the stream
% eps: iid packet loss rate
% =====================================================================

k = k(:); in = in(:);
N = length(k); % number of frames in the intra-period

% First column: frame arrival probabilities
% Second column: frame arrival probabilities with an additional FEC block
phi = zeros(N,2);

if templayers > 1
    mode = 1; % hP
else
    mode = 0; % IPPP
end

for i = 1:N
    phi(i,1) = sum(binopdf( 0:in(i)  , in(i)+k(i)  , eps));
    phi(i,2) = sum(binopdf( 0:in(i)+1, in(i)+k(i)+1, eps));
end

if f == 0
    best_score = calc_score( phi(:,1), mode );
    best_alloc = zeros(N,1);
    return
end

best_alloc = in;
best_score = calc_score( phi(:,1), mode );

for i = 1:f % for each additional FEC block
    
    % check each neighbor
    for j = 1:N
        
        phi_    = phi(:,1);
        phi_(j) = phi(j,2);
        
        score = calc_score( phi_, mode );
        
        if score > best_score
            best_score = score;
            picked = j;
        end
    end
    
    best_alloc(picked) = best_alloc(picked)+1;
    
    % need to only change the success prob corr. to the picked frame
    phi(picked,1) = phi(picked,2);
    
    % to update the value in the second column
    phi(picked,2) = sum( binopdf(0:(best_alloc(picked)+1), k(picked)+best_alloc(picked)+1, eps) );
end

return