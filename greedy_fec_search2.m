function [best_alloc, best_score] = greedy_fec_search2( f, k, in, templayers, varargin )
% GREEDY_FEC_SEARCH2    greedy search for distributing f FEC packets on
% frames of length k (should be faster than GREEDY_FEC_SEARCH3 for independent losses)
% f: total number of FEC blocks we can use
% k: vector of frame lengths in the intra-period 
% in: initial FEC block distribution
% alpha: transition probability from good to bad
% beta: transition probability from bad to good
% templayers: number of temporal layers of the stream
% =====================================================================

k = k(:); in = in(:);
N = length(k); % number of frames in the intra-period
independent = 0;

% First column: frame arrival probabilities
% Second column: frame arrival probabilities with an additional FEC block
phi = zeros(N,2);

if templayers > 1
    mode = 1; % hP
else
    mode = 0; % IPPP
end

if length(varargin) == 1
    
    independent = 1;
    eps = varargin{1};
    
    for i = 1:N
        phi(i,1) = sum(binopdf( 0:in(i)  , in(i)+k(i)  , eps));
        phi(i,2) = sum(binopdf( 0:in(i)+1, in(i)+k(i)+1, eps));
    end   
else
    
    alpha = varargin{1};
    beta  = varargin{2};
    n = max(k)+f+1;
    Pi = alpha/(alpha+beta);
    
    R = create_mtx( n, [1,beta*(1-alpha).^(0:n-2)], [1-beta,alpha*beta*(1-alpha).^(0:n-2)] );
    S = create_mtx( n, [1,alpha*(1-beta).^(0:n-2)], [1-alpha,beta*alpha*(1-beta).^(0:n-2)] );
    
    for i = 1:N
            M = in(i); K = k(i);
        if M > 0
            phi(i,1) = (1-Pi) * sum(S( K:K+M  , K+M   )) + Pi * sum(R( 1:M  , K+M   ));
        else
            phi(i,1) = (1-Pi) * sum(S( K:K+M  , K+M   ));
        end
        if f > 0
            phi(i,2) = (1-Pi) * sum(S( K:K+M+1, K+M+1 )) + Pi * sum(R( 1:M+1, K+M+1 ));
        end
    end
end

if f == 0
    best_score = calc_score( phi(:,1), mode );
    best_alloc = zeros(N,1);
    return
end

% best_alloc = zeros(N,1); % initial best allocation: no FEC blocks on frames
best_alloc = in; % initial best allocation: no FEC blocks on frames
% best_score = 0;          % initial best score = 0
best_score = calc_score( phi(:,1), mode );

for i = 1:f % for each additional FEC block
    
    % check each neighbor
    for j = 1:N
        
        phi_ = phi(:,1); phi_(j) = phi(j,2);
        
        score = calc_score( phi_, mode );
        
        if score > best_score
            best_score = score;
            picked = j;
        end
    end
    
    best_alloc(picked) = best_alloc(picked)+1;
    
    % need to only change the success prob corr. to the picked frame
    phi(picked,1) = phi(picked,2);
    
    M = best_alloc(picked)+1; % to update the value in the second column
    K = k(picked);
    if independent    
        phi(picked,2) = sum( binopdf(0:M, K+M, eps) );
    else
        phi(picked,2) = (1-Pi) * sum(S( K:K+M, K+M )) + Pi * sum(R( 1:M, K+M ));
    end
end

return

function A = create_mtx( n, B, b )
% Make sure B and b have correct length (row vectors of size(A,1))
% B: row vector st B(i)=Pr(0^(i-1) | 1) or Pr(1^(i-1) | 0)
% b: row vector st b(i)=Pr(0^(i-1)1 | 1) or Pr(1^(i-1)0 | 0)
A = zeros(n);
A(1,:) = B;
for i = 2:n
    for j = i:n
        A(i,j) = fliplr(b(1:(j-i+1))) * A(i-1,i-1:j-1)';
    end
end

return