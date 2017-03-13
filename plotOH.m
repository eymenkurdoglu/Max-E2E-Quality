%% Initialize
close all

f0 = figure; hold all; box on;

sequences = {'CREW', 'CITY', 'HARBOUR', 'SOCCER'};

%% Plot
for j = 1 : length(sequences)
    
    video = sequences{j};
    
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
                
        if L == 1
            QSTAR_IPP = max( mnqq( vs.q0 .* ( (bw./vs.R0).*(15/vs.fmax)^-vs.beta_f ).^(-1/vs.beta_q), vs.alpha_q, vs.qmin )...
    .* mnqt( 15, vs.alpha_f, vs.fmax ), mnqq( vs.q0 .* ( (bw./vs.R0) ).^(-1/vs.beta_q), vs.alpha_q, vs.qmin ));
            plot( bw/1e3, QSTAR_IPP, [color,line], 'LineWidth', 2 );
        else
            QSTAR_hPP = max( mnqq( vs.q0 .* ( (bw./vs.R0).*(15/vs.fmax)^-vs.beta_f ).^(-1/vs.beta_q), vs.alpha_q, vs.qmin )...
    .* mnqt( 15, vs.alpha_f, vs.fmax ), mnqq( vs.q0 .* ( (bw./vs.R0) ).^(-1/vs.beta_q), vs.alpha_q, vs.qmin ));
            plot( bw/1e3, QSTAR_hPP, [color,line], 'LineWidth', 2 );
            QSTAR = (QSTAR_IPP-QSTAR_hPP)./QSTAR_IPP;
        end
    end
end

ylim( [0.3 1.01] )
xlim( [0.1 1.6] )
xlabel('Encoding Bitrate (Mbps)');
ylabel('Encoded Video Quality')
legend('CREW-IPP','CREW-hPP','CITY-IPP','CITY-hPP','HARBOUR-IPP','HARBOUR-hPP',...
    'SOCCER-IPP','SOCCER-hPP','Location','SouthEast');

saveTightFigure(f0,'~/Google Drive/NYU/Research/papers/TMM-2/fig/codingOH.eps')
clear all