%% Initialize
close all

f0 = figure; 
% f1 = figure; hold all; box on;

sequences = {'CREW', 'CITY', 'HARBOUR', 'SOCCER'};
% sequences = {'CREW'};

%% Plot
for j = 1 : length(sequences)
    
    video = sequences{j};
    subplot(1,2,1); hold all; box on;
    for L = [1 3]
        resfile = ['results/',video,'-',num2str(L),'.mat'];
        if exist( resfile, 'file' ); load( resfile )
        else fprintf( 'Result file for %s with %d layers missing!\n', video, L ); return; end
        
        if      strcmp(video,'CREW'); color = 'r';
        elseif  strcmp(video,'CITY'); color = 'g';
        elseif  strcmp(video,'SOCCER'); color = 'b';
        elseif  strcmp(video,'HARBOUR'); color = 'k';
        elseif  strcmp(video,'ICE'); color = 'm';
        elseif  strcmp(video,'FG'); color = 'c';
        end
        
        if L == 1; line = '-'; else line = '--'; end
    
        q_range_15 = vs.q0 .* ( (bw./vs.R0).*(15/vs.fmax)^-vs.beta_f ).^(-1/vs.beta_q);
        q_range_30 = vs.q0 .* ( (bw./vs.R0) ).^(-1/vs.beta_q);
        NQQ_15 = mnqq( q_range_15, vs.alpha_q, vs.qmin );
        NQQ_30 = mnqq( q_range_30, vs.alpha_q, vs.qmin );
        NQT_15 = mnqt( 15, vs.alpha_f, vs.fmax );

        Q_range = 0.1 : 0.01 : 1;
        intres = 1 - (Q_range*(1-exp(-vs.alpha_q))/NQT_15);
        intres( intres<0 ) = min( intres( intres>0 ) );
        q_15 = -(vs.alpha_q*vs.qmin)./log( intres );
        q_30 = -(vs.alpha_q*vs.qmin)./log( 1-Q_range*(1-exp(-vs.alpha_q)) );        
        
        if L == 1
            QSTAR_IPP = max( NQQ_15.* NQT_15, NQQ_30 );

            RATE_IPP = vs.R0*min( 2^-vs.beta_f .* (vs.q0./q_15).^vs.beta_q, (vs.q0./q_30).^vs.beta_q );

            plot( bw/1e3, QSTAR_IPP, [color,line], 'LineWidth', 2 );
        else            
            QSTAR_hPP = max( NQQ_15.* NQT_15, NQQ_30 );

            RATE_hPP = vs.R0*min( 2^-vs.beta_f .* (vs.q0./q_15).^vs.beta_q, (vs.q0./q_30).^vs.beta_q );            
            
            plot( bw/1e3, QSTAR_hPP, [color,line], 'LineWidth', 2 );
        end
    end
    subplot(1,2,2); hold all; box on;
    plot(0.1:0.01:1, 100*((RATE_hPP./RATE_IPP)-1), [color,'.'], 'LineWidth', 2)
end
% 100*((RATE_hPP./RATE_IPP)-1)
subplot(1,2,1)
ylim( [0.3 1.01] )
xlim( [0.1 1.6] )
xlabel('Encoding Bitrate (Mbps)');
ylabel('Encoded Video Quality')
legend('CREW-IPP','CREW-hPP','CITY-IPP','CITY-hPP','HARBOUR-IPP','HARBOUR-hPP',...
    'SOCCER-IPP','SOCCER-hPP','Location','SouthEast');

subplot(1,2,2)
xlim( [0.1 1] )
ylabel('Coding Overhead Percentage');
xlabel('Encoded Video Quality')

% saveTightFigure(f0,'~/Google Drive/NYU/Research/papers/TMM-2/fig/codingOH.eps')
clear all