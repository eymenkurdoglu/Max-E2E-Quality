clc
close all
dbstop if error
path = '~/Dropbox/Matlab/UEP4hP/data/';

%% data process
for folder = dir(path)'
    if folder.isdir && ~( strcmp(folder.name,'.') || strcmp(folder.name,'..') )
        if ~exist( [folder.name,'.mat'], 'file' )
            
            fprintf( 'Could not find prior simulations, simulating for %s\n', folder.name );
            packet_size = 200; % bytes
            alphas = 10.^(-1:-1:-3);
            betas = 9*alphas;
            I = 4;
            
            param = strsplit(folder.name,'-');
            sequence   = param{1};
            spatialres = param{2};
            fps        = str2double( param{3} );
            N          = str2double( param{4} ); range = ((I-1)*N + 1) : I*N;  T = N/fps;
            templayers = str2double( param{5} );
            l = []; qp = []; rates = [];
            for file = dir( strcat(path,folder.name,'/*.txt') )'
                contents = dlmread( strcat(path,folder.name,'/',file.name) );
                l  = [l contents(:,1)];
                qp = [qp contents(:,2)];
                dum = strsplit(file.name,'.');
                rates = [rates str2double(dum{1})];
            end
            [rates,dum] = sort(rates);
            qp = qp(range,dum);
            l = ceil( l(range,dum) / packet_size );
            clear range contents file dum rates
            abw = max(sum(l));
            fec = abw-sum(l);
            src_perc = 100*sum(l)/abw;

            [fec_allocs, max_qual, max_fps] = simulate( alphas, betas, fec, l, qp, T, templayers );
            save( [folder.name,'.mat'] )
        else
            load( [folder.name,'.mat'] )
        end
    end
end

%% plot and learn
figure; hold all;
for j = 1:length(alphas)
    plot( src_perc, max_fps(j,:) );
end
ylabel('Expected received fps')
xlabel('Source rate percentage (%)')

figure; hold all;
for j = 1:length(alphas)
    plot( src_perc, max_qual(j,:) );
end
ylabel('Quality (MNQQ*MNQT)')
xlabel('Source rate percentage (%)')

figure;
title('FEC rate per frame')
opt_fec_rate = [];
for j = 1:length(alphas)
    Q = max_qual(j,:);
    opt_fec_alloc = squeeze( fec_allocs(:,Q==max(Q),j) );
    opt_src_lengt = l(:,Q==max(Q));
    subplot(3,1,j);
    bar(0:N-1, opt_fec_alloc./(opt_fec_alloc+opt_src_lengt))
    ylim( [0 1] )
    xlim( [-1 N] )
    if j==1
        title('FEC rate per frame')
    elseif j ==3
        xlabel('Frame index')
    end
    legend(['P_{GB} = ',num2str(alphas(j))],'Location','NorthEast')
end
