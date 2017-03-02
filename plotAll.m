function redundancyRatios = plotAll( sequences, L )

willPlotNQQandNQT = 0;
willPlotMonteCarlo = 0;
willSavePlotsAsEps = 0;
isMarkovian = 1;

channel = '';
if isMarkovian; channel = '-Markov'; end

%% 2x2 plots
close all
f1 = figure; f4 = figure; f5 = figure; f6 = figure;
if willPlotNQQandNQT;  f2 = figure; f3 = figure; end
if willPlotMonteCarlo; f7 = figure; f8 = figure; f9 = figure; f10 = figure; end
f11 = figure;

for j = 1 : length(sequences)
    
    video = sequences{j};
    
    load(['results/',video,'-',num2str(L),channel,'.mat'])
    
    epsilon = pgb./(pgb+pbg);
    epsilon = epsilon(epsilon <= 0.2);

    numChains = length(epsilon);
    numCapacs = length(bw);
    
    lambda = 1./pbg(1:numChains);
    
    xl = [0 bw(numCapacs)]/1e3;
    
    LEGEND = cell(1,numChains); % legend strings
    for i = 1:numChains
        if isMarkovian
            LEGEND{i} = ['\lambda = ',num2str(lambda(i))];
        else
            LEGEND{i} = ['\epsilon = ',num2str(pgb(i))];
        end
    end
    
%% ##### Q-STAR #####
    figure(f1)
        subplot(2,2,j)
        hold all
        box on
        xlabel( 'Sending Bitrate (Mbps)' )
        title( video )
        xlim( xl )
%         ylim( [0.2 1] )
        
        plot( bw/1e3, (NQQ(1:numChains,:).*NQT(1:numChains,:))' )
        
        if ~isMarkovian
            
            Q_lossless = max( mnqq( vs.q0 .* ( (bw./vs.R0).*(15/vs.fmax)^-vs.beta_f ).^(-1/vs.beta_q)...
                , vs.alpha_q, vs.qmin ).* mnqt( 15, vs.alpha_f, vs.fmax ),mnqq( vs.q0 .* ( (bw./vs.R0)...
                ).^(-1/vs.beta_q), vs.alpha_q, vs.qmin ) );
            
            plot( bw/1e3, Q_lossless );  
            
            if j == 3; legend( [LEGEND,'\epsilon = 0'], 'Location', 'best' ); end
        else
            if j == 3; legend( LEGEND, 'Location', 'best' ); end
        end
    
    if willPlotNQQandNQT
        figure(f2); subplot(2,2,j);
        plot(bw/1e3,NQQ(1:numChains,:)'); ylabel('NQQ');
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(LEGEND,'Location',...
                'Best'); end

        figure(f3); subplot(2,2,j);
        plot(bw/1e3,NQT(1:numChains,:)'); ylabel('NQT');
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(LEGEND,'Location',...
                'Best'); end 
    end
    
%% ##### TOTAL FEC REDUNDANCY #####
    
    TotalFECPerc = 100 * (1-(R(1:numChains,:)./repmat(bw,numChains,1))');
    
    figure(f4)
        subplot(2,2,j)
        hold all
        box on
        xlabel('Sending Bitrate (Mbps)')
        title(video)
        xlim( xl )
        ylim( [0 100] )
        
        plot( bw/1e3, TotalFECPerc )
        if j == 3; legend(LEGEND,'Location','Best'); end

%% ##### OPTIMIZED FRAME RATES #####
    figure(f5)
        subplot(2,2,j)
        hold all
        box on
        xlabel( 'Sending Bitrate (Mbps)' )
        title( video )
        xlim( xl )
        ylim( [0 31] ); 
        
        plot( bw/1e3, F(1:numChains,:)' )
    
        if ~isMarkovian
            
            FR_lossless = 15+15*(mnqq( vs.q0 .* ( (bw./vs.R0).*(15/vs.fmax)^-vs.beta_f )...
                .^(-1/vs.beta_q), vs.alpha_q, vs.qmin ).* mnqt( 15, vs.alpha_f, vs.fmax ) < mnqq( vs.q0 ...
                .* ( (bw./vs.R0) ).^(-1/vs.beta_q), vs.alpha_q, vs.qmin ));
            
            plot( bw/1e3, FR_lossless )
            
            if j == 3; legend( [LEGEND,'\epsilon = 0'], 'Location', 'Best' ); end
        else
            if j == 3; legend( LEGEND, 'Location', 'Best' ); end
        end
    
%% ##### FEC REDUNDANCY RATIOS PER LAYER #####

    redundancyRatios = zeros(numChains,numCapacs,L+1);
    for v = 1:numChains
        for u = 1:numCapacs
            k = K{v,u};
            m = M{v,u};
            redundancyRatios(v,u,:) = 100*find_mean_per_layer( m./(k+m), length(m), L );
        end
    end

%     avgRedundancyRatios = squeeze( mean(redundancyRatios,2) );
    figure(f11)
        for v = 1:numChains
            subplot( length(sequences), numChains, (j-1)*numChains+v )
            plot( bw/1e3, squeeze( redundancyRatios(v,:,:) ) )
            xlim( xl )
%             ylim( [0 55] )
            if j == 1; title( ['$\epsilon$ = ',num2str(epsilon(v))], 'Interpreter', 'latex' ); end
            if v == 1; ylabel( video, 'FontSize', 12 ); end
        end
    
    figure(f6)
        subplot(2,2,j)
        hold all
        box on
        title(video)
        
        affineModels = zeros(2,1);
        horAxis = 100*epsilon;

        model = fit( horAxis', mean(TotalFECPerc)', 'poly1' );
        x = linspace(min(horAxis),max(horAxis),200);
        affineModels(1) = plot( x, x * model.p1 + model.p2 ); % Overall (fit)
        scatter( horAxis, mean(TotalFECPerc) ); % Overall (actual)
        
        load(['results/',video,'-1.mat'])
        TotalFECPerc = 100 * (1-(R(1:numChains,:)./repmat(bw,numChains,1))');
        model = fit( horAxis', mean(TotalFECPerc)', 'poly1' );
        x = linspace(min(horAxis),max(horAxis),200);
        affineModels(2) = plot( x, x * model.p1 + model.p2 ); % Overall (fit)
        scatter( horAxis, mean(TotalFECPerc) ); % Overall (actual)        
            
        xlabel( 'Packet Loss Rate \epsilon (%)' )
        if j == 3; legend( affineModels, 'hPP', 'IPP', 'Location', 'Best' ); end
        
    if willPlotMonteCarlo && exist([video,'-',num2str(L),'-MC.mat'],'file')
        load([video,'-',num2str(L),'-MC.mat'])

        figure(f7); subplot(2,2,j);
        plot(bw/1e3,30*MeanFrameIntervals(1:numChains,:)'); ylabel('Mean Frame Interval'); ylim([0 15])
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(LEGEND,'Location',...
                'Best'); end

        figure(f8); subplot(2,2,j);
        plot(bw/1e3,30*StdFrameIntervals(1:numChains,:)'); ylabel('Std Frame Interval'); ylim([0 10])
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(LEGEND,'Location',...
                'Best'); end

        figure(f9); subplot(2,2,j);
        plot(bw/1e3,numFreezes(1:numChains,:)'/100e3); ylabel('Pr(Frz)'); ylim([0 0.025])
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(LEGEND,'Location',...
                'Best'); end

        figure(f10); subplot(2,2,j);
        plot(bw/1e3,MeanNumDecFrames(1:numChains,:)'); ylabel('MNDF');
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(LEGEND,'Location',...
                'Best'); end
    end
    
end

if willSavePlotsAsEps
    target = '~/Google Drive/NYU/Research/papers/fec/fig/';
    saveTightFigure(f1,[target,'quality-',num2str(L),channel,'.eps'])
    saveTightFigure(f4,[target,'videoBitrate-',num2str(L),channel,'.eps'])
    saveTightFigure(f5,[target,'encFrRate-',num2str(L),channel,'.eps'])
    saveTightFigure(f6,[target,'affine.eps'])
%     saveTightFigure(f6,[target,'fecRatesPerLayer-',num2str(L),channel,'.eps'])
end

return

function mean_x = find_mean_per_layer( x, N, L )
% find_mean_per_layer       groups the x vector according to the
% intra-period structure
%  This function creates cells for each layer + I frames. Then, the entries
%  of the x vector are put inside the corresponding cells. 

mean_x = cell(1,L+1); % put frame sizes per each layer in their corresponding cell

% P-frames from enhancement layers
for l = L : -1 : 2
    mean_x{ l+1 } = x(2 : 2 : end);
    x(2 : 2 : end) = [];
end

% I-frames
mean_x{1} = x(1 : N*2^(1-L) : end);
x(1 : N*2^(1-L) : end) = [];

% P-frames from base layer
mean_x{2} = x;

mean_x = (cellfun(@mean, mean_x))';

return