function k = estimFrameSz( vs, R )
% estimFrameSz      estimate the sizes of the I and P frames in the
% intra-period given the target encoding bitrate
%  Size of the P-frame tau ms apart from its reference normalized wrt the
%  I-frame is estimated using the model: | zhat | = 1 - exp(-theta*tau^-eta)
%  

    f = vs.f;
    L = vs.L;
    N = f*vs.ipr; % number of frames in intra-period
    B = 1000 * vs.ipr * R/8; % byte budget in intra-period
    
    eta   = vs.eta  ( vs.fr == vs.f );
    theta = vs.theta( vs.fr == vs.f );
    
    zhatNorm = fliplr( 1-exp(-theta*(1000*(2.^(0:L-1))./f).^eta) );

    % n: number of P-frames in intra-period per TL, e.g. [7 8 16]
    n = fliplr( N ./ (2.^[ 1:L-1, L-1]) );
    n(1) = n(1)-1;

    zhat = (B(:)/(1 + n*zhatNorm')*[1, zhatNorm])';

    % find frame indices for each TL and create the k vector or matrix
    % start from the last layer, then previous, ...
    k = zeros(N,length(B));
    for r = 1:length(B)
        for layer = L : -1 : 1   
            k(mod(0:N-1,2^(L-layer))==0,r) = zhat(layer+1,r);
        end
    end
    
    k(1,:) = zhat(1,:);
    
return