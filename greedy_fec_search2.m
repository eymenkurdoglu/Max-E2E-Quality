function [best_alloc, best_score] = greedy_fec_search2( f, k, varargin )
% Faster greedy FEC search. This function finds a distribution of "f" FEC
% blocks over the frames of size k, such that the expected number of
% decodable frames is maximized. [best_alloc, best_score] =
% greedy_fec_search2( f, k, p ) f: total number of FEC blocks we can use k:
% vector of frame sizes in the intra-period eps: packet loss rate (assuming
% uncorrelated losses)

k = k(:);
N = length(k); % number of frames in the intra-period

mode = 1; % IPPP

if length(varargin) == 1
    independent = 1;
    eps = varargin{1};
elseif length(varargin) == 2
    independent = 0;
    p = varargin{1};
    q = varargin{2};
end

% this table (Nx2) holds the success probabilities for each frame. 
% First column: success probs with the current distribution of FEC blocs. 
% Second column: success probs with 1 additional FEC block for each frame
if independent % independent losses
    success = [ binopdf(0,k,eps), binopdf(0,k+1,eps)+binopdf(1,k+1,eps) ];
elseif f > 0 % markovian losses
    success = zeros(N,2);
    n = max(k)+f;
    R = create_mtx( n, [1,q*(1-p).^(0:n-2)], [1-q,q*p*(1-p).^(0:n-2)] );
    S = create_mtx( n, [1,p*(1-q).^(0:n-2)], [1-p,q*p*(1-q).^(0:n-2)] );
    for i = 1:N
        success(i,1) = q*S(k(i),k(i))/(p+q);
        success(i,2) = ( q*sum(S(k(i):k(i)+1,k(i)+1)) + p*R(1,k(i)+1) )/(p+q);
    end
    success(:,1)
else
    success = zeros(N,1);
    n = max(k);
    S = create_mtx( n, [1,p*(1-q).^(0:n-2)], [1-p,q*p*(1-q).^(0:n-2)] );
    for i = 1:N
        success(i) = q*S(k(i),k(i))/(p+q);
    end
    best_score = calc_score( success, mode );
    best_alloc = zeros(N,1);
    return
end

best_alloc = zeros(N,1); % initial best allocation: no FEC blocks on frames
best_score = 0; % initial best score = 0

for i = 1:f % for each additional FEC block
    
    for j = 1:N % try each frame in the intra-period
        success_ = success(:,1);
        success_(j) = success(j,2);
        score = calc_score( success_, mode );
        if score > best_score
            best_score = score;
            picked = j;
        end
    end
    
    best_alloc(picked) = best_alloc(picked)+1;
    % need to only change the success prob corr. to the picked frame
    success(picked,1) = success(picked,2);
    
    M = best_alloc(picked)+1; % to update the value in the second column
    K = k(picked);
    if independent    
        success(picked,2) = sum( binopdf(0:M, K+M, eps) );
    else
        success(picked,2) = q/(p+q)*sum(S(K:K+M, K+M)) + p/(p+q)*sum(R(1:M,K+M));
    end
end

return