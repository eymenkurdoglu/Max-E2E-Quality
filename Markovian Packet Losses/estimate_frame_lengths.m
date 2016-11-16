function k = estimate_frame_lengths(FR, L, N, B, video )
% FR: encoding frame rate, L: number of temporal layers, N: number of
% frames in intra-period, B: byte budget in intra-period

if      strcmp(video, 'CREW')
% CREW, FR=30, a=0.093498, b=0.433992
% CREW, FR=15, a=0.161280, b=0.423905
    if FR == 15 
        a = 0.161280; b = 0.423905;
    elseif FR == 30
        a = 0.093498; b = 0.433992;
    end
elseif  strcmp(video, 'CITY')
% CITY, FR=30, a=0.121674, b=0.304659
% CITY, FR=15, a=0.183831, b=0.322710
    if FR == 15 
        a = 0.183831; b = 0.322710;
    elseif FR == 30
        a = 0.121674; b = 0.304659;
    end
elseif  strcmp(video, 'HARBOUR')
% HARBOUR, FR=30, a=0.088230, b=0.405442
% HARBOUR, FR=15, a=0.166287, b=0.370140
    if FR == 15 
        a = 0.166287; b = 0.370140;
    elseif FR == 30
        a = 0.088230; b = 0.405442;
    end
elseif  strcmp(video, 'FG')
% FG, FR=30, a=0.039767, b=0.576068
% FG, FR=15, a=0.060339, b=0.402718
    if FR == 15
        a = 0.060339; b = 0.402718;
    elseif FR == 30
        a = 0.039767; b = 0.576068;
    end    
end

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