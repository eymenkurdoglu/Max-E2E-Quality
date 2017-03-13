close all
L = 3;
channel = '-Markov';
sequences = {'CREW', 'CITY', 'HARBOUR', 'SOCCER'};
figure;
b = bw==1600;
for j = 1 : length(sequences)
    video = sequences{j};
    load(['results/',video,'-',num2str(L),channel,'.mat'])
    epsilon = pgb./(pgb+pbg);
    epsilon = epsilon(epsilon <= 0.2);
    numChains = length(epsilon);
    numCapacs = length(bw);
    lambda = 1./pbg(1:numChains);
    redundancyRatios = zeros(numChains,L+1);
    for v = 1:numChains
        k = K{v,b};
        m = M{v,b};
        redundancyRatios(v,:) = 100*find_mean_per_layer( m./(k+m), length(m), L );
    end
    subplot(2,2,j)
    hold all
    box on
    xlabel( 'Mean burst duration' )
    title( video )   
    F(:,b)
    plot( 1.6 .* lambda' .* F(:,b)./1600, redundancyRatios ) ;
%     xlim( [0, max(B)] );
    ylim( [0 100] );
end