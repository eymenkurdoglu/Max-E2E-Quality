close all
clear all
L = 3;
channel = '-LRvLmbd';
sequences = {'CITY', 'HARBOUR', 'CREW', 'SOCCER'};
f = figure;
for j = 1 : length(sequences)
    video = sequences{j};
    load( ['results/15Hz/',video,'-',num2str(L),channel,'.mat'] )
    epsilon = pgb./(pgb+pbg);
    numChains = length(epsilon);
    numCapacs = length(bw);
    lambda = 1./pbg(1:numChains);
    redundancyRatios = zeros(numChains,L+1);
    for v = 1:numChains
        k = K{v};
        m = M{v};
        redundancyRatios(v,:) = 100*find_mean_per_layer( m./(k+m), length(m), L );
    end
    if j <= 2
        subplot(2,2,j+2)
        hold all
        box on
        xlabel( 'Mean burst duration (\lambda)' )
        title( video )
        plot( lambda, redundancyRatios' ) ;
        xlim( [min(lambda), max(lambda)] );
        ylim( [0 100] );
        if j == 2; legend('I-frame','TL(1)','TL(2)','TL(3)','Location','SouthWest'); end
    end
end
saveTightFigure( '/Users/eymen/Google Drive/NYU/Research/papers/TMM-2/fig_/layerRedunvsLambda.eps' )