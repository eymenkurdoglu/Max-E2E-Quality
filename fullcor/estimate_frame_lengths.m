function k = estimate_frame_lengths(FR, L, N, B )
% FR: encoding frame rate, L: number of temporal layers, N: number of
% frames in intra-period, B: byte budget in intra-period

switch FR
    case 15
        a = 0.1613;
    case 30
        a = 0.0935;
    otherwise
        display('Unable to estimate frame lengths, terminating.');
        return
end
b = 0.429;

% rho = P-frame length / I-frame length (assumed to be constant regardless of the bitrate)
rho = fliplr( 1-exp(-1*a*(1000*(2.^(0:L-1))./FR).^b) );

% h: how many frames in intra-period per layer, eg. [7 8 16]
h = fliplr(  N ./ (2.^[ 1:L-1, L-1]) );
h(1) = h(1) - 1;

% frame lengths in bytes for each TL at the given bitrate
len = (B(:)/(1 + h*rho')*[1, rho])';

% find frame indices for each TL and create the k vector or matrix
% start from the last layer, then previous, ...
k = zeros(N,length(B));
for r = 1:length(B)
    for layer = L : -1 : 1   
        k(mod(0:N-1,2^(L-layer))==0,r) = len(layer+1,r);
    end
end
k(1,:) = len(1,:);
return