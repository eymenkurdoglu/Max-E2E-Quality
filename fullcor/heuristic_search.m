function [Qmax, topt, Ropt, mopt] = heuristic_search( PLR, R_S, L, video )
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

%% RSTAR and QSTAR parameters for each video
if      strcmp(video, 'CREW')
    beta_q=1.061165; beta_t=0.707289; Rmax=1188.956800; qmin=14.100638;
    alpha_q = 4.51; alpha_t = 3.09; 
elseif  strcmp(video, 'CITY')
    beta_q=1.142239; beta_t=0.471515; Rmax=1150.621600; qmin=10.662216;
    alpha_q = 7.25; alpha_t = 4.10; 
elseif  strcmp(video, 'HARBOUR')
    beta_q=1.320098; beta_t=0.584034; Rmax=1196.192000; qmin=21.765459;
    alpha_q = 9.65; alpha_t = 2.83; 
elseif  strcmp(video, 'ICE')    
    beta_q=0.868087; beta_t=0.653480; Rmax=943.812800; qmin=5.318390;
    alpha_q = 5.61; alpha_t = 3.00; 
elseif  strcmp(video, 'FOREMAN')    
    beta_q=1.087917; beta_t=0.558064; Rmax=1171.495200; qmin=9.878532;
    alpha_q = 4.57; alpha_t = 3.80;
end
tmax = 30;

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
        fprintf('Trying video rate=%d, M=%d\n',R,Total - sum(k) - sum(m));
        
        % now find the maximum mean NQT we can achieve with M FEC packets
        % and the allocation that achieves this
        % WILL CHANGE
        [m, NQT] = greedy_fec_search( Total - sum(k) - sum(m), k, tree, PLR, alpha_t, tmax, RBR_IPR, m );
        
        % estimated mean QS at this bitrate
        q = qmin * ( (RBR_VIDEO_BITRATES(j)/Rmax)*(t/tmax)^(-1*beta_t) )^(-1/beta_q);
        
        % calculate the subjective perceptual quality at this STAR
        Q = mnqq(q,alpha_q,qmin) * NQT;
        
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