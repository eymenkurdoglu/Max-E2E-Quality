path = '~/Dropbox/Matlab/3_OptimalQuality/data/';
% videos = {};
% handle = figure;

%% data process
for folder = dir(path)'
    if folder.isdir && ~( strcmp(folder.name,'.') || strcmp(folder.name,'..') )
        if ~exist( ['./',folder.name,'.mat'], 'file' )
            
            fprintf( 'Could not find prior simulations, simulating for %s\n', folder.name );
            packet_size = 200; % bytes
            alphas = [0.1, 0.05, 0.02, 0.01, 0.005, 0.001];
            betas = 9*alphas;
            I = 4;
            
            param      = strsplit(folder.name,'-');
            vid        = param{1}; if strcmp( vid, 'AKIYO' ) continue; end
            spatialres = param{2};
            fps        = str2double( param{3} ); if fps ~= 30 continue; end
            N          = str2double( param{4} ); range = ((I-1)*N + 1) : I*N;  T = N/fps;
            templayers = str2double( param{5} );
            
            lengths = []; qp = []; rates = [];
            for file = dir( strcat(path,folder.name,'/*.txt') )'
                contents = dlmread( strcat(path,folder.name,'/',file.name) );
                lengths  = [lengths contents(:,1)];
                qp = [qp contents(:,2)];
                dum = strsplit(file.name,'.');
                rates = [rates str2double(dum{1})];
            end
            [rates,dum] = sort(rates);
            qp = qp(range,dum);
            lengths = ceil( lengths(range,dum) / packet_size );
            clear range contents file dum rates param
            abw = max(sum(lengths));
            fec = abw-sum(lengths);
            src_perc = 100*sum(lengths)/abw;

            [fec_allocs, max_qual, max_fps] = simulate( alphas, betas, fec, lengths, qp, T, templayers, vid );
            
            save( [folder.name,'.mat'] )
        else
            load( [folder.name,'.mat'] )
            if fps ~= 30 continue; end
            videos = [videos, vid];
            optimized_qual = zeros(size(alphas));
            opt_src_perc = zeros(size(alphas));
            opt_qp = zeros(size(alphas));
            opt_fps = zeros(size(alphas));
            for j = 1:length(alphas)
                ix = max_qual(j,:) == max(max_qual(j,:));
                opt_src_perc(j) = src_perc( ix );
                optimized_qual(j) = max(max_qual(j,:));
                opt_qp(j) = mean(qp(:,ix));
                opt_fps(j) = max_fps(j,ix);
            end
            
            figure; hold all;
            for j = 1:length(alphas)
                plot( src_perc, max_fps(j,:) );
            end
            title(vid)
            ylabel('Expected received fps')
            xlabel('Source rate percentage (%)')

            figure; hold all;
            for j = 1:length(alphas)
                plot( src_perc, max_qual(j,:) );
            end
            title(vid)
            ylabel('Quality (MNQQ*MNQT)')
            xlabel('Source rate percentage (%)')
            
%             subplot(2,2,1); hold all; box on;
%             plot( 1./betas, opt_src_perc ); %ylim([30 90]);
%             ylabel('SRP* (%)'); xlabel('Avg. burst length');
%             subplot(2,2,2); hold all; box on;
%             plot( 1./betas, optimized_qual ); %ylim([0.6 1]);
%             ylabel('Qual*'); xlabel('Avg. burst length');
%             subplot(2,2,3); hold all; box on;
%             plot( 1./betas, opt_qp ); 
%             ylabel('QP*'); xlabel('Avg. burst length'); %ylim([15 40])
%             subplot(2,2,4); hold all; box on;
%             plot( 1./betas, opt_fps );
%             xlabel('Avg. burst length'); ylabel('FPS*'); %ylim([25 31])
        end
    end
end
legend(videos,'Location','Best')
saveTightFigure('all_seq.eps')
%% plot and learn

% figure;
% title('FEC rate per frame')
% opt_fec_rate = [];
% for j = 1:length(alphas)
%     Q = max_qual(j,:);
%     opt_fec_alloc = squeeze( fec_allocs(:,Q==max(Q),j) );
%     opt_src_lengt = l(:,Q==max(Q));
%     subplot(3,1,j);
%     bar(0:N-1, opt_fec_alloc./(opt_fec_alloc+opt_src_lengt))
%     ylim( [0 1] )
%     xlim( [-1 N] )
%     if j==1
%         title('FEC rate per frame')
%     elseif j ==3
%         xlabel('Frame index')
%     end
%     legend(['P_{GB} = ',num2str(alphas(j))],'Location','NorthEast')
% end
