function [best_fec, maxqual, max_num_frames] = heuristic_search( alpha, beta, fec, lengths, q, templayers )

global alpha_q minQS alpha_t maxt ipr

N = size(lengths,1); % intra-period length
maxqual = 0;
max_num_frames = 0;
best_fec = [];

ix = fec >= 0;

if ~any(ix)
    return
end

lengths = lengths(:,ix);
q = q(ix);
fec = fec(ix);

best_distribution_so_far = zeros(N,1);

for i = length(fec):-1:1
    
    num_fec_blocks = fec(i)-sum(best_distribution_so_far);

    if alpha+beta == 1
        [best_distribution_so_far, max_num_frames] = greedy_fec_search2( num_fec_blocks, lengths(:,i), best_distribution_so_far, templayers, alpha );
    else
        [best_distribution_so_far, max_num_frames] = greedy_fec_search3( fec(i), lengths(:,i), alpha, beta, templayers );
    end

    Q = mnqq(q(i),alpha_q,minQS) * mnqt(max_num_frames/ipr,alpha_t,maxt);

    if Q > maxqual 
        maxqual = Q;
    else
%         fprintf('Max quality found\n')
        best_fec = i+1; % the one before, so i+1
        return
    end
end
best_fec = i;
return

function Q = mnqq( q, alpha_q, qmin )

num = 1-exp(-1*alpha_q*(qmin./q));
denom = 1-exp(-1*alpha_q);

Q = num/denom;

return

function Q = mnqt( t, alpha_t, tmax )

num = 1-exp(-1*alpha_t*(t./tmax).^0.63);
denom = 1-exp(-1*alpha_t);

Q = num/denom;

return