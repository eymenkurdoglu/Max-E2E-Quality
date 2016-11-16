function [Qmax, topt, Ropt, mopt, PMFopt] = heuristic( PLR, R_S, L, video )

RBR_IPR = 16/15;
RBR_PACKET_SIZE = 200;
RBR_FRAME_RATES = [30, 15];

RBR_VIDEO_BITRATES = R_S : -1 : 75;
RBR_SEARCH_PERSISTENCE = 4;

% max number of packets we can send at this sending rate
Total = floor( R_S*RBR_IPR / (0.008*RBR_PACKET_SIZE) );

%% Start the search
Qmax = 0;

for t = RBR_FRAME_RATES

    N = t * RBR_IPR;
    
    Param = video_param( video, RBR_IPR, L, N ); tree = create_hierP_tree( L, N );
    
    m = zeros(N,1); highest_Q_so_far = 0; best_pmf_so_far = 0;
    
    best_rate_so_far = R_S*(1-2*PLR);
    k = estimate_frame_lengths(t, L, N, RBR_IPR*best_rate_so_far/0.008, video);
    k = ceil( k / RBR_PACKET_SIZE );
    M = Total - sum(k);
    [m, NQT, pmf] = heuristic_fec( M, k, tree, PLR, Param );
    q = Param.qmin * ( (R/Param.Rmax)*(t/Param.tmax)^-Param.beta_t )^(-1/Param.beta_q);
    NQQ = mnqq(q,Param.alpha_q,Param.qmin);
        Q = NQQ * NQT;
        fprintf('Q = %4.4f x %4.4f = %4.4f\n',NQQ,NQT,Q);
        
    if highest_Q_so_far > Qmax
        Qmax = highest_Q_so_far;
        PMFopt = best_pmf_so_far;
        topt = t;
        Ropt = best_rate_so_far;
        mopt = best_alloc_so_far;
    end
end

return