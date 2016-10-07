function [Qmax, topt, Ropt, mopt] = heuristic_search( alpha, beta, R_S, L, video )
% Inputs: PLR: iid packet loss rate, R_S: Sending rate, FR: frame rate
% Outputs: Qmax: max achievable perceptual quality, topt: maximizer frame
% rate, Ropt: maximizer video bitrate, mopt: maximizer FEC strategy

RBR_IPR = 16/15;
RBR_PACKET_SIZE = 200;
RBR_FRAME_RATES = [30, 15];
RBR_VIDEO_BITRATES = R_S : -1 : 50;
RBR_SEARCH_PERSISTENCE = 4;

% max number of packets we can send at this sending rate
Total = floor( R_S*RBR_IPR / (0.008*RBR_PACKET_SIZE) );

param = star_param( video, RBR_IPR );

%% Start the search
Qmax = 0;
topt = 0;
Ropt = 0;

for i = 1:length(RBR_FRAME_RATES)
    
    t = RBR_FRAME_RATES(i);
    N = RBR_IPR * t;
    tree = create_hierP_tree( L, N );
    
    % estimated mean frame lengths (in packets) in each TL at this bitrate
    FrameSizes = ceil( estimate_frame_lengths(t, L, RBR_IPR, RBR_VIDEO_BITRATES) / RBR_PACKET_SIZE );
    
    prevk = zeros(N,1);
    m = zeros(N,1);
    highest_Q_so_far = 0;
    best_rate_so_far = 0;
    descending = 0;
    
    for j = 1:length(RBR_VIDEO_BITRATES)
        
        k = FrameSizes(:,j);
        if all( k == prevk ) || Total < sum(k)
            continue; 
        end
        prevk = k;
        fprintf('Trying video rate=%d, M=%d\n',RBR_VIDEO_BITRATES(j),Total - sum(k) - sum(m));
        
        M = Total - sum(k) - sum(m);
        markov = create_prob_struct( alpha, beta, max(k) + M + 1 );
        % now find the maximum mean NQT we can achieve with M FEC packets
        % and the allocation that achieves this
        [m, NQT] = greedy_fec_search3( M, k, tree, markov, param, m );
        
        % estimated mean QS at this bitrate
        q = param.qmin * ( (RBR_VIDEO_BITRATES(j)/param.Rmax)*(t/param.tmax)^(-1*param.beta_t) )^(-1/param.beta_q);
        
        % calculate the subjective perceptual quality at this STAR
        Q = mnqq(q,param.alpha_q,param.qmin) * NQT;
        
        if Q > highest_Q_so_far
            descending = 0;
            highest_Q_so_far = Q;
            best_rate_so_far = RBR_VIDEO_BITRATES(j);
            best_alloc_so_far = m;
        else
            descending = descending + 1;
            if descending > RBR_SEARCH_PERSISTENCE
                break;
            end
        end
    end
    
    if highest_Q_so_far > Qmax
        Qmax = highest_Q_so_far;
        topt = t;
        Ropt = best_rate_so_far;
        mopt = best_alloc_so_far;
    end
end

return