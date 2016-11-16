% script to show how t_d/t_e changes with respect to video bitrate, FEC
% bitrate, and packet loss rate

videos = {'CREW'};

for i = 1 : length(videos)
    
    for PLR = 5:5:30
        
        path = [videos{i},'/PLR_',num2str(PLR),'%.mat'];

        if exist( path, 'file' )

            load( path )
        else
            RBR_FRAME_RATES = [30,15]; %Hz
            bitrateRange = 25:25:1250; %kbps
            RBR_NUM_TEMP_LYRS = 3;
            RBR_PACKET_SIZE = 200; %bytes
            RBR_IPR = 16/15; %sec

            t_d = zeros(length(RBR_FRAME_RATES),length(bitrateRange),length(bitrateRange));

            for j = 1:length( RBR_FRAME_RATES )

                FR = RBR_FRAME_RATES(j);
                N = RBR_IPR * FR;

                % update this
                FrameSizes = frame_lengths(FR, RBR_NUM_TEMP_LYRS, RBR_IPR, bitrateRange, RBR_PACKET_SIZE);

                for u = 1:length(bitrateRange)

                    k = FrameSizes(:,u);

                    for v = 1:length( bitrateRange )

                        M = ceil( bitrateRange(v) * RBR_IPR / (0.008*RBR_PACKET_SIZE) );

                        [~, score] = greedy_fec_search( M, k, zeros(N,1), RBR_NUM_TEMP_LYRS, PLR/100 );

                        t_d( j, u, v ) = score/RBR_IPR;

                        if abs((t_d( j, u, v )/FR)-1) <= 1e-4
                            t_d( j, u, v+1 : end ) = FR;
                            break;
                        end
                    end
                end
            end
            save( path )
        end
        
        figure;
        
        subplot(1,length(RBR_FRAME_RATES),1)
        imagesc(flipud(squeeze(t_d(1,:,:))));
        xlabel('FEC rate'); ylabel('Video rate'); title('FR=30 Hz')
        colorbar

        subplot(1,length(RBR_FRAME_RATES),2)
        imagesc(flipud(squeeze(t_d(2,:,:))));
        xlabel('FEC rate'); ylabel('Video rate'); title('FR=15 Hz')
        colorbar
    end
end