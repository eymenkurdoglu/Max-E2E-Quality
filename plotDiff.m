function plotDiff( eps, markov, montecarlo, varargin )
dbstop if error
%% Initialize
close all

separateqstar = 0;
if nargin > 3; separateqstar = varargin{4}; end

if      markov == 0;    channel = '';
elseif  markov == 1;    channel = '-Markov';
else    display Wrong data identifier!; return;
end

f1 = figure; f2 = figure; f9 = figure;
if separateqstar;   f3 = figure; f4 = figure; end
if montecarlo;      f5 = figure; f6 = figure; f7 = figure; f8 = figure; end

sequences = {'CREW', 'CITY', 'HARBOUR', 'SOCCER'};

for j = 1 : length(sequences)
    
    video = sequences{j};
    
    affineModel = zeros(2,1);
    
    %% Prepare
    for L = [1 3]
        resfile = ['results/',video,'-',num2str(L),channel,'.mat'];
        if      exist( resfile, 'file' ); load( resfile )
        else    fprintf( 'Result file for %s with %d layers missing%s!\n', video, L, channel ); return; end
    
        epsilon = pgb./(pgb+pbg);
        epsilon = epsilon(epsilon <= 0.2); % greedy search ineffective out of this range
        lambda = 1./pbg(epsilon <= 0.2);
        numChains = length(epsilon);
        numCapacs = length(bw);
    
        l = cell(1,numChains); % legend strings
        if      markov == 1; for i = 1:numChains; l{i}=['\lambda = ',num2str(lambda(i))];   end
        elseif  markov == 0; for i = 1:numChains; l{i}=['\epsilon = ',num2str(epsilon(i))]; end
        end
        
        D = D(1:numChains,1:numCapacs); F = F(1:numChains,1:numCapacs);
        K = K(1:numChains,1:numCapacs); M = M(1:numChains,1:numCapacs); R = R(1:numChains,1:numCapacs);
        NQQ = NQQ(1:numChains,1:numCapacs); NQT = NQT(1:numChains,1:numCapacs);        
        
        if L == 1; diffQ   = NQQ' .* NQT';  else diffQ = (diffQ - NQQ' .* NQT')./diffQ; end
        if L == 1; diffNQQ = NQQ';          else diffNQQ = diffNQQ - NQQ';              end
        if L == 1; diffNQT = NQT';          else diffNQT = diffNQT - NQT';              end
        
        % Plot affine models
        FECPerc = 100 * (1 - R ./ repmat(bw,numChains,1));
        model = fit( 100*epsilon', mean(FECPerc,2), 'poly1' );
        x = 100 * linspace( min(epsilon), max(epsilon), 200 );
        
        figure( f9 ); subplot(2,2,j); hold all; box on
        affineModel( (L+1)/2 ) = plot( x, x * model.p1 + model.p2 ); % Overall (fit)
        scatter( 100*epsilon', mean(FECPerc,2) ) % Overall (actual)
        fprintf([video,'-',num2str(L),': %fx+%f\n'],model.p1,model.p2);
        
        if L == 1; diffFECPerc = FECPerc; else diffFECPerc = diffFECPerc - FECPerc; end
        
        if montecarlo
            mcresfile = ['results/',video,'-',num2str(L),channel,'-MC.mat'];
            if exist( mcresfile, 'file' ); load( mcresfile )
            else fprintf( 'MonteCarlo result file for %s with %d layers missing%s!\n', video, L, channel ); return; end
            
            MeanFrameIntervals = MeanFrameIntervals(1:numChains,:); StdFrameIntervals = StdFrameIntervals(1:numChains,:);
            numFreezes = numFreezes(1:numChains,:); MeanNumDecFrames = MeanNumDecFrames(1:numChains,:);

            if L == 1; MFI = MeanFrameIntervals';   else MFI = (MFI - MeanFrameIntervals'); end
            if L == 1; SFI = StdFrameIntervals';    else SFI = (SFI - StdFrameIntervals');  end
            if L == 1; FRP = numFreezes'/1e5;       else FRP = (FRP - numFreezes'/1e5);     end
            if L == 1; MNDF = MeanNumDecFrames';    else MNDF = (MNDF - MeanNumDecFrames'); end
        end
    end
    
    xlabel( 'Packet Loss Rate \epsilon (%)' )
    legend( affineModel, 'hPP', 'IPP', 'Location', 'Best' );
    
    xl = [bw(1) bw(numCapacs)]/1e3;
    
    %% Plot others
    figure(f1); subplot(2,2,j); hold all; box on;
    plot( bw/1e3, 100*diffQ );
    xlabel( 'Sending Bitrate (Mbps)' ); title( video ); xlim( xl );
    if j == 1; legend(l,'Location','NorthEast'); end
    
    figure(f2); subplot(2,2,j);
    plot( bw/1e3, diffFECPerc' );
    xlabel( 'Sending Bitrate (Mbps)' ); title( video ); xlim( xl );
    if j == 3; legend(l,'Location','NorthEast'); end
    
    if separateqstar   
        figure(f3); subplot(2,2,j);
        plot( bw/1e3, diffNQQ ); 
        xlabel( 'Sending Bitrate (Mbps)' ); title( video ); xlim( xl );
        if j == 3; legend(l,'Location','NorthEast'); end

        figure(f4); subplot(2,2,j);
        plot( bw/1e3, diffNQT );
        xlabel( 'Sending Bitrate (Mbps)' ); title( video ); xlim( xl );
        if j == 3; legend( l, 'Location', 'SouthEast' ); end
    end
    
    if montecarlo
        figure(f5); subplot(2,2,j); hold all; box on;
        for v = 1:numChains
            scatter( bw/1e3, 1000*MFI(:,v), '*' );
        end
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); ylim( [-30 200] )
        if j == 3; legend(l,'Location','NorthEast'); end

        figure(f6); subplot(2,2,j); hold all; box on;
        for v = 1:numChains
            scatter( bw/1e3, 1000*SFI(:,v), '*' );
        end
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); ylim( [-30 200] )
        if j == 3; legend(l,'Location','NorthEast'); end

        figure(f7); subplot(2,2,j); hold all; box on;
        for v = 1:numChains
            scatter( bw/1e3, FRP(:,v), '*' );
        end
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl);
        if j == 3; legend(l,'Location','NorthEast'); end

        figure(f8); subplot(2,2,j); hold all; box on;
        for v = 1:numChains
            scatter( bw/1e3, MNDF(:,v), '*' );
        end        
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); 
        if j == 3; legend( l, 'Location', 'SouthEast' ); end
    end
end
%% Save
if eps
    target = '~/Google Drive/NYU/Research/papers/TMM-2/fig/';
    saveTightFigure(f1,[target,'quality-diff',channel,'.eps'])
    saveTightFigure(f2,[target,'videoBitrate-diff',channel,'.eps'])
    if montecarlo
        saveTightFigure(f5,[target,'avgFrInt-diff.eps'])
        saveTightFigure(f6,[target,'stdFrInt-diff.eps'])
    end
end

return