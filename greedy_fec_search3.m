function [best_alloc, best_score] = greedy_fec_search3( f, k, alpha, beta, templayers )
% Faster greedy FEC search. This function finds a distribution of "f" FEC
% blocks over the frames of size k, such that the expected number of
% decodable frames is maximized. [best_alloc, best_score] =
% greedy_fec_search2( f, k, p ) f: total number of FEC blocks we can use k:
% vector of frame sizes in the intra-period eps: packet loss rate (assuming
% uncorrelated losses)

k = k(:);
N = length(k); % number of frames in the intra-period
n = max(k)+f+1;

mode = 0;
if templayers > 1
    mode = 1;
end

% "0" and "1" duration tail distribution functions and probability mass functions
P_R = [1,beta*(1-alpha).^(0:n-2)];
P_L = [1,alpha*(1-beta).^(0:n-2)];
p_R = [1-beta,beta*alpha*(1-alpha).^(0:n-2)];
p_L = [1-alpha,beta*alpha*(1-beta).^(0:n-2)];

L  = create_mtx( n, P_R, p_R ); % Pr(m-1 losses occur in the next n-1 following a loss)
L0 = create_mtx( n, p_R, p_R ); % Pr(m-1 losses occur in the next n-1 between 2 losses)
L1 = L-L0; % Pr(m-1 losses occur in the next n-1 between a loss and a reception)

R  = create_mtx( n, P_L, p_L ); % Pr(m-1 receptions occur in the next n-1 following a reception)
R1 = create_mtx( n, p_L, p_L ); % Pr(m-1 receptions occur in the next n-1 between 2 receptions)
R0 = R-R1; % Pr(m-1 receptions occur in the next n-1 between a reception and a loss)

T = [1-beta, alpha; beta, 1-alpha];

best_alloc = zeros(N,1); % initial best allocation: no FEC blocks on frames
pr = framearrprobs( L0, L1, R0, R1, T, k, best_alloc );
best_score = calc_score( pr, mode );

for i = 1:f % for each additional FEC block
    
    curr_alloc = best_alloc;
    
    for j = 1:N
        
        neighbor = curr_alloc;
        neighbor(j) = neighbor(j) + 1;
        
        score = calc_score( framearrprobs( L0, L1, R0, R1, T, k, neighbor ), mode );
        
        if score > best_score
            best_score = score;
            best_alloc = neighbor;
        end
        
    end
    
end

return