function [best_alloc, best_score] = greedy_fec_search3( f, k, alpha, beta, templayers )
% GREEDY_FEC_SEARCH3    greedy search for distributing f FEC packets on
% frames of length k
% f: total number of FEC blocks we can use
% k: vector of frame lengths in the intra-period 
% alpha: transition probability from good to bad
% beta: transition probability from bad to good
% templayers: number of temporal layers of the stream
% =====================================================================
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

T = [1-beta,   alpha;
       beta, 1-alpha];

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

function Phi = framearrprobs( L0, L1, R0, R1, T, k, m )
% 
N = length(k);
Phi = zeros(2,N);

for i = 1:N
    
    % find the reference frame, ref<0 iff i==1
    z = mod(i,4);
    if z == 1
        ref = i-4;
    elseif z == 2 || z == 0
        ref = i-1;
    else
        ref = i-2;
    end
    
    A = [ sum( L0(1 : m(i)-1, k(i)+m(i)-1) ), sum( R0(k(i) : k(i)+m(i)-1, k(i)+m(i)-1) );
          sum( L1(1 : m(i), k(i)+m(i)-1) ), sum( R1(k(i)-1 : k(i)+m(i)-1, k(i)+m(i)-1) ) ];
      
    if ref < 0;
        Phi(:,1) = A * [T(1,2);T(2,1)]/sum(T(1,2)+T(2,1));
    else
        inb = ref:i; inb(1) = []; inb(end) = [];
        Phi(:,i) = A * T^(sum(k(inb)+m(inb))+1) * Phi(:,ref)/sum(Phi(:,ref));
    end
    
end
Phi = (sum(Phi))';
return