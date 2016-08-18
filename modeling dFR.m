clc
% script to show how t_d/t_e changes wrt. R_v, R_fec, and PLR
RBR_VIDEO = 'CREW';
RBR_FRAME_RATES = {'30.0','15.0'}; %Hz
RBR_TARGET_BITRATES = 25:25:1250; %kbps
RBR_NUM_TEMP_LYRS = 3;
RBR_PACKET_SIZE = 200; %bytes
RBR_IPR = 16/15; %sec
RBR_SELECTED_IPR = 4;

plr = 0.3;

l = zeros( 1, RBR_NUM_TEMP_LYRS+1 );
a = [ 0.0935 0.1613 ]; b = 0.429;

t_d = zeros(length(RBR_FRAME_RATES),length( RBR_TARGET_BITRATES ),length( RBR_TARGET_BITRATES ));

for fr = 1:length( RBR_FRAME_RATES )
    FR = str2double( RBR_FRAME_RATES{fr} );
    N = RBR_IPR * FR;
    k = zeros(N,1);
    G = 2^(RBR_NUM_TEMP_LYRS-1);
    n = fliplr(  N ./ (2.^[ 1:RBR_NUM_TEMP_LYRS-1, RBR_NUM_TEMP_LYRS-1]) ); % number of frames in IPR per layer
    n(1) = n(1) - 1; % account for the I-frame
    temp_distance = 1000*(2.^(0:RBR_NUM_TEMP_LYRS-1))./FR; %msec
    rho = fliplr( 1-exp(-1*a(fr)*temp_distance.^b) ); % P-frame length / I-frame length
    for rv = 1:length( RBR_TARGET_BITRATES )
        l(1) = RBR_TARGET_BITRATES(rv) * RBR_IPR / (1 + n * rho');
        l(2:end) = l(1) * rho; l = l/0.008;
        k(1:G:end) = l(2); k(1) = l(1); k(3:G:end) = l(3); k(2:2:end) = l(4);
        k = ceil( k / RBR_PACKET_SIZE );
        for rf = 1:length( RBR_TARGET_BITRATES )
            M = ceil( RBR_TARGET_BITRATES(rf) * RBR_IPR / (0.008*RBR_PACKET_SIZE) );
            fprintf('Optimizing, %d FEC packets, %d kbps\n',M,RBR_TARGET_BITRATES(rv));
            [~, score] = greedy_fec_search( M, k, zeros(N,1), RBR_NUM_TEMP_LYRS, plr );
            t_d( fr, rv, rf ) = score/RBR_IPR;
            if abs((t_d( fr, rv, rf )/FR)-1) <= 1e-4
                fprintf('!');
                t_d( fr, rv, rf+1 : end ) = FR;
                break;
            end
        end
    end
end
save('t_d_030.mat')

figure;
imagesc(flipud(squeeze(t_d(1,:,:))));
xlabel('FEC rate'); ylabel('Video rate'); title('FR=30 Hz')
colorbar
figure;
imagesc(flipud(squeeze(t_d(2,:,:))));
xlabel('FEC rate'); ylabel('Video rate'); title('FR=15 Hz')
colorbar