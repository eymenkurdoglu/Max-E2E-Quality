clc
% script to show how t_d/t_e changes wrt. R_v, R_fec, and PLR
RBR_VIDEO = 'CREW';
RBR_FRAME_RATES = {'30.0','15.0'}; %Hz
RBR_TARGET_BITRATES = 25:25:1250; %kbps
RBR_NUM_TEMP_LYRS = 3;
RBR_PACKET_SIZE = 200; %bytes
RBR_IPR = 16/15; %sec

RBR_PLR = '030';
plr = str2double( RBR_PLR )/100;

if ~exist( ['plot_t_d/CREW/t_d_',RBR_PLR,'.mat'], 'file' )

    t_d = zeros(length(RBR_FRAME_RATES),length( RBR_TARGET_BITRATES ),length( RBR_TARGET_BITRATES ));

    for fr = 1:length( RBR_FRAME_RATES )

        FR = str2double( RBR_FRAME_RATES{fr} );
        N = RBR_IPR * FR;

        FrameSizes = frame_lengths(FR, RBR_NUM_TEMP_LYRS, RBR_IPR, RBR_VIDEO_BITRATES, RBR_PACKET_SIZE);

        for rv = 1:length( RBR_TARGET_BITRATES )

            k = FrameSizes(:,rv);

            for rf = 1:length( RBR_TARGET_BITRATES )

                M = ceil( RBR_TARGET_BITRATES(rf) * RBR_IPR / (0.008*RBR_PACKET_SIZE) );

                [~, score] = greedy_fec_search( M, k, zeros(N,1), RBR_NUM_TEMP_LYRS, plr );

                t_d( fr, rv, rf ) = score/RBR_IPR;

                if abs((t_d( fr, rv, rf )/FR)-1) <= 1e-4
                    t_d( fr, rv, rf+1 : end ) = FR;
                    break;
                end
            end
        end
    end
    save(['plot_t_d/CREW/t_d_',RBR_PLR,'.mat'])
else
    load(['plot_t_d/CREW/t_d_',RBR_PLR,'.mat'])
end

figure;
imagesc(flipud(squeeze(t_d(1,:,:))));
xlabel('FEC rate'); ylabel('Video rate'); title('FR=30 Hz')
colorbar

figure;
imagesc(flipud(squeeze(t_d(2,:,:))));
xlabel('FEC rate'); ylabel('Video rate'); title('FR=15 Hz')
colorbar