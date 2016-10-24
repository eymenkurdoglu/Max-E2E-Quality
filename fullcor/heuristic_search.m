function [Qmax, topt, Ropt, mopt, PMFopt] = heuristic_search( p_bg, p_gb, R_S, L, video )
% INPUTS: p_bg, p_gb: 2-state loss mechanism parameters, R_S: Sending rate,
% L: number of temporal layers, video: name of the video sequence, OUTPUTS:
% Qmax: max achievable perceptual quality, topt: maximizer frame rate,
% Ropt: maximizer video bitrate, mopt: maximizer FEC strategy

RBR_IPR = 16/15;
RBR_PACKET_SIZE = 200;
RBR_FRAME_RATES = [30, 15];
RBR_VIDEO_BITRATES = R_S : -1 : 75;
RBR_SEARCH_PERSISTENCE = 4;

% max number of packets we can send at this sending rate
Total = floor( R_S*RBR_IPR / (0.008*RBR_PACKET_SIZE) );

%% Start the search
Qmax = 0; topt = 0; Ropt = 0;

for t = RBR_FRAME_RATES
    fprintf('# Trying encoding frame rate = %d fps...\n',t);
    N = t * RBR_IPR;
    Param = video_param( video, RBR_IPR, L, N );
    
    % estimated mean frame lengths in packets at this frame rate and bitrate
    k = ceil( estimate_frame_lengths(t, L, N, ...
        RBR_IPR*RBR_VIDEO_BITRATES/0.008) / RBR_PACKET_SIZE );
    
    % init
    Markov = create_prob_struct( p_bg, p_gb, Total );
    prevk = zeros(N,1); descending = 0;
    m = zeros(N,1); highest_Q_so_far = 0; best_rate_so_far = 0; best_pmf_so_far = 0;

    for i = 1:length(RBR_VIDEO_BITRATES)
        
        R = RBR_VIDEO_BITRATES(i);
        if all( k(:,i) == prevk ) || Total < sum(k(:,i))
            continue; 
        end
        prevk = k(:,i);
        
        M = Total - sum(k(:,i)) - sum(m);
        fprintf('## Trying video rate = %d kbps, M = %d ==> ',R,M);       
        
        % maximize mean NQT with M FEC packets
        tic; [m, NQT, pmf] = greedy_fec_search( M, k(:,i), Markov, Param, m ); toc;
        
        % estimated mean QS at this bitrate
        q = Param.qmin * ( (R/Param.Rmax)*(t/Param.tmax)^-Param.beta_t )^(-1/Param.beta_q);
        
        % calculate the subjective perceptual quality at this STAR
        Q = mnqq(q,Param.alpha_q,Param.qmin) * NQT;
        
        if Q > highest_Q_so_far
            descending = 0;
            best_pmf_so_far = pmf;
            highest_Q_so_far = Q;
            best_rate_so_far = R;
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
        PMFopt = best_pmf_so_far;
        topt = t;
        Ropt = best_rate_so_far;
        mopt = best_alloc_so_far;
    end
end

return