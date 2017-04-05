function layerRedun = plotAll( eps, L, markov, montecarlo, varargin )
%% Initialize
close all

separateqstar = 0;
oneFR = '';
if nargin > 4; separateqstar = varargin{1}; end
if nargin > 5; oneFR = [num2str(varargin{2}),'Hz/']; end

if      markov == 0;    channel = '';
elseif  markov == 1;    channel = '-Markov';
else    display Wrong data identifier!; return;
end

f1 = figure; f2 = figure; f3 = figure; f4 = figure; f11 = figure;
if separateqstar;   f5 = figure; f6 = figure; end
if montecarlo;      f7 = figure; f8 = figure; f9 = figure; f10 = figure; end

sequences = {'CITY', 'HARBOUR', 'CREW', 'SOCCER'};

for j = 1 : length(sequences)
    
    video = sequences{j};
    
    resfile = ['results/',oneFR,video,'-',num2str(L),channel,'.mat'];
    if      exist( resfile, 'file' ); load( resfile )
    else    fprintf( 'Result file for %s with %d layers missing%s!\n', video, L, channel ); return; end
    
    epsilon = pgb./(pgb+pbg);
    epsilon = epsilon(epsilon <= 0.2); % greedy search ineffective out of this range
    lambda = 1./pbg(epsilon <= 0.2);
    numChains = length(epsilon);
    numCapacs = length(bw);
    xl = [bw(1) bw(numCapacs)]/1e3;
    
    D = D(1:numChains,1:numCapacs); F = F(1:numChains,1:numCapacs);
    K = K(1:numChains,1:numCapacs); M = M(1:numChains,1:numCapacs); R = R(1:numChains,1:numCapacs);
    NQQ = NQQ(1:numChains,1:numCapacs); NQT = NQT(1:numChains,1:numCapacs);
    
    l = cell(1,numChains); % legend strings
    if      markov == 1; for i = 1:numChains; l{i}=['\lambda = ',num2str(lambda(i))];   end
    elseif  markov == 0; for i = 1:numChains; l{i}=['\epsilon = ',num2str(epsilon(i))]; end
    end
    
    TotalFECPercentage = 100 * ( 1 - (R./repmat(bw,numChains,1)) );
    
%% ##### Q-STAR #####
    figure(f1); subplot(2,2,j); hold all; box on
    plot( bw/1e3, NQQ' .* NQT' )
    xlabel( 'Sending Bitrate (Mbps)' ); title( video ); xlim( xl );
    
    if markov == 0
        Q_lossless = max( mnqq( vs.q0 .* ( (bw./vs.R0).*(15/vs.fmax)^-vs.beta_f ).^(-1/vs.beta_q)...
            , vs.alpha_q, vs.qmin ).* mnqt( 15, vs.alpha_f, vs.fmax ),mnqq( vs.q0 .* ( (bw./vs.R0)...
            ).^(-1/vs.beta_q), vs.alpha_q, vs.qmin ) );
        plot( bw/1e3, Q_lossless );  
        if j == 4; legend( [l,'\epsilon = 0'], 'Location', 'best' ); end
    else
        if j == 4; legend( l, 'Location', 'best' ); end
    end
    
    if separateqstar
        figure(f5); subplot(2,2,j);
        plot( bw/1e3, NQQ' );
        xlabel('Sending Bitrate (Mbps)'); ylabel('NQQ'); title(video); xlim(xl); 
        if j == 4; legend(l,'Location','Best'); end

        figure(f6); subplot(2,2,j);
        plot( bw/1e3, NQT' );
        xlabel( 'Sending Bitrate (Mbps)' ); ylabel( 'NQT' ); title( video ); xlim( xl );
        if j == 4; legend(l,'Location','Best'); end
    end    
%% ##### TOTAL FEC REDUNDANCY #####   
    figure(f2); 
    if j <= 2
        subplot(2,2,j+2); hold all; box on
        if markov == 1; plot( bw/24, TotalFECPercentage' ); xlabel('Avg. Block Size (packets)' ); xlim( 1000*xl/24 ); end
        if markov == 0; plot( bw/1e3, TotalFECPercentage' ); xlabel('Sending Bitrate (Mbps)' ); xlim( xl ); end
        title(video); ylim( [0 100] ); 
        if j == 2; legend(l,'Location','Best'); end
    end
%% ##### OPTIMIZED FRAME RATES #####
    figure(f3); subplot(2,2,j); hold all; box on
    plot( bw/1e3, F(1:numChains,:)' )
    xlabel( 'Sending Bitrate (Mbps)' ); title( video ); xlim( xl ); ylim( [0 31] ); 

    if markov == 0
        FR_lossless = 15+15*(mnqq( vs.q0 .* ( (bw./vs.R0).*(15/vs.fmax)^-vs.beta_f )...
            .^(-1/vs.beta_q), vs.alpha_q, vs.qmin ).* mnqt( 15, vs.alpha_f, vs.fmax ) < mnqq( vs.q0 ...
            .* ( (bw./vs.R0) ).^(-1/vs.beta_q), vs.alpha_q, vs.qmin ));
        plot( bw/1e3, FR_lossless )
        if j == 3; legend( [l,'\epsilon = 0'], 'Location', 'Best' ); end
    else
        if j == 3; legend( l, 'Location', 'Best' ); end
    end
%% ##### AVERAGE dFR #####
    dFR = zeros(numChains,numCapacs);
    for v = 1:numChains
        for u = 1:numCapacs
            dFR(v,u) = (0:(vs.ipr*F(v,u)))*D{v,u}'/vs.ipr;
        end
    end
    figure(f11); subplot(2,2,j); hold all; box on
    plot( bw/1e3, dFR' )
    xlabel( 'Sending Bitrate (Mbps)' ); title( video ); xlim( xl ); ylim( [0 31] );     
%% ##### FEC REDUNDANCY RATIOS PER LAYER #####
    layerRedun = zeros(numChains,numCapacs,L+1);
    for v = 1:numChains
        for u = 1:numCapacs; k = K{v,u}; m = M{v,u};
            layerRedun(v,u,:) = 100*find_mean_per_layer( m./(k+m), length(m), L );
        end
    end

    figure(f4)
    for v = 1:numChains
        if j <= 2
            subplot( 2, numChains, (j-1)*numChains+v )
            plot( bw/1e3, squeeze( layerRedun(v,:,:) ) )
            xlim( xl ); ylim( [0 100] )
            if v == 1; ylabel( video, 'FontSize', 12 ); end
            if markov == 1; if j == 1; title( ['$\lambda$ = ',num2str(lambda(v))], 'Interpreter', 'latex' ); end
            else if j == 1; title( ['$\epsilon$ = ',num2str(epsilon(v))], 'Interpreter', 'latex' ); end
            end
        end
    end
        
    if montecarlo
        mcresfile = ['results/',oneFR,video,'-',num2str(L),channel,'-MC.mat'];
        if exist( mcresfile, 'file' ); load( mcresfile )
        else fprintf( 'MonteCarlo result file for %s with %d layers missing%s!\n', video, L, channel ); return; end
        
        MeanFrameIntervals = MeanFrameIntervals(1:numChains,:); StdFrameIntervals = StdFrameIntervals(1:numChains,:);
        numFreezes = numFreezes(1:numChains,:); MeanNumDecFrames = MeanNumDecFrames(1:numChains,:);        

        figure(f7); subplot(2,2,j); hold all; box on;
        for v = 1:numChains
            scatter( bw/1e3, 1000*MeanFrameIntervals(v,:), '*' );
        end        
        ylabel('Mean Frame Interval');
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl);
        if j == 4; legend(l,'Location','Best'); end

        figure(f8); subplot(2,2,j); hold all; box on;
        for v = 1:numChains
            scatter( bw/1e3, 1000*StdFrameIntervals(v,:), '*' );
        end            
        ylabel('Std Frame Interval');
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl);
        if j == 4; legend(l,'Location','Best'); end
        
        figure(f9); subplot(2,2,j); hold all; box on
        for v = 1:numChains
            scatter( bw/1e3, numFreezes(v,:)/1e5, '*' );
        end
        ylabel('Pr(Frz)'); xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl);
        if j == 4; legend(l,'Location','Best'); end

        figure(f10); subplot(2,2,j); hold all; box on;
        for v = 1:numChains
            scatter( bw/1e3, MeanNumDecFrames(v,:), '*' );
        end         
        ylabel('MNDF'); xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl);
        if j == 4; legend(l,'Location','Best'); end
    end
    
end
%% Save
if eps
    target = '~/Google Drive/NYU/Research/papers/TMM-2/fig_/';
    saveTightFigure(f1,[target,'quality-',num2str(L),channel,'.eps'])
    saveTightFigure(f2,[target,'videoBitrate-',num2str(L),channel,'.eps'])
%     saveTightFigure(f3,[target,'encFrRate-',num2str(L),channel,'.eps'])
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