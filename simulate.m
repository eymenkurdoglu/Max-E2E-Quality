function [fec_allocs, max_qual, max_fps] = simulate( alphas, betas, fec, lengths, qp, T, templayers, vid )
% SIMULATE  simulate Gilbert loss scenarios for an intra-period in a video stream.
% alphas    : vector of from-good-to-bad probabilites
% betas     : vector of from-bad-to-good probabilities
% fec       : vector of FEC packet amounts
% lengths   : frame-length column vectors for each FEC packet amount in fec
% qp        : QP-value column vectors for each FEC packet amount in fec
% T         : intra-period duration in seconds
% templayers: number of temporal layers in the stream (works only for =3 or =1 for now)
%======================================================================================

N = size( lengths,1 ); % intra-period length
fec_allocs = zeros(N,length(fec),length(alphas));
max_qual   = zeros(length(alphas),length(fec));
max_fps    = zeros(length(alphas),length(fec));

% The parameter ?t controls how fast the NQT drops as t decreases
% with a smaller value corresponding to a faster dropping rate.
% The parameter ?q controls how fast the NQQ drops as q increases
% with a smaller value corresponding to a slower dropping rate.
qmin = 10; tmax = 30;
if      strcmp(vid, 'CREW')
    alpha_q = 4.51; alpha_t = 3.09; 
elseif  strcmp(vid, 'CITY')
    alpha_q = 7.25; alpha_t = 4.10; 
elseif  strcmp(vid, 'HARBOUR')
    alpha_q = 9.65; alpha_t = 2.83; 
elseif  strcmp(vid, 'ICE')    
    alpha_q = 5.61; alpha_t = 3.00; 
elseif  strcmp(vid, 'FOREMAN')    
    alpha_q = 4.57; alpha_t = 3.80;
end

for j = 1:length(alphas)

    alpha = alphas(j);
    beta  = betas(j);
    avg_num_frames  = zeros(1,length(fec));
    
    fprintf('\nSimulating for good-to-bad prob = %f\n', alpha);

    % search!
    for i = length(fec):-1:1
        tic
%         if alpha+beta == 1
            [best_alloc, max_num_frames] = greedy_fec_search2( fec(i), lengths(:,i), zeros(N,1), templayers, alpha, beta );
%         else
%             [best_alloc, max_num_frames] = greedy_fec_search3( fec(i), lengths(:,i), alpha, beta, templayers );
%         end
        toc
        avg_num_frames(i) = max_num_frames;
        fec_allocs(:,i,j) = best_alloc;
    end
    max_fps(j,:) = avg_num_frames/T;
    max_qual(j,:) = mnqq(mean(qp),alpha_q,qmin) .* mnqt(avg_num_frames/T,alpha_t,tmax);
end

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