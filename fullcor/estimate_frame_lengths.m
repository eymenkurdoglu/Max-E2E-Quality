function k = estimate_frame_lengths(FR, L, RBR_IPR, rv )

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

N = FR * RBR_IPR;

rv = rv(:);
S = length(rv);

% rho = P-frame length / I-frame length (assumed to be constant regardless of the bitrate)
rho = fliplr( 1-exp(-1*a*(1000*(2.^(0:L-1))./FR).^b) );

% n holds the numbers of frames in one IPR for each layer, eg. [7 8 16]
n = fliplr(  N ./ (2.^[ 1:L-1, L-1]) );
n(1) = n(1) - 1;

% frame lengths in bytes for each TL at the given bitrate
l = (rv*RBR_IPR/0.008 / (1 + n * rho'))*[1, rho];

k = zeros(N,S);

for t = 1:S
    i = 0;
    for j = fliplr( l(t,2:end) )
        k(mod(0:N-1,2^i)==0,t) = j;
        i = i+1;
    end

    k(1,t) = l(t,1);
end

return