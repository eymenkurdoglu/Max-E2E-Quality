function plotAll( sequences, L )

willPlotNQQandNQT = 0;
willPlotMonteCarlo = 0;
willSavePlotsAsEps = 1;
thisIsMarkovian = 0;
fourSeq = 0;

channel = '';
if thisIsMarkovian; channel = '-Markov'; end
if length(sequences) > 1; fourSeq = 1; end

%% 2x2 plots, section a
close all
f1 = figure; f4 = figure; f5 = figure; f6 = figure;
if willPlotNQQandNQT;  f2 = figure; f3 = figure; end
if willPlotMonteCarlo; f7 = figure; f8 = figure; f9 = figure; f10 = figure; end

for j = 1 : length(sequences)
    
    video = sequences{j};
    
    load([video,'-',num2str(L),channel,'.mat'])
    
    epsilon = pgb./(pgb+pbg);
    epsilon = epsilon(epsilon <= 0.2);

    numChains = length(epsilon);
    numCapacs = length(bw);
    
    lambda = 1./pbg(1:numChains);
    
    xl = [0 bw(numCapacs)]/1e3;
    
    TotalFECPerc = 100*(1-(R(1:numChains,:)./repmat(bw,numChains,1))');
    
    l = cell(1,numChains); % legend strings
    for i = 1:numChains
        if thisIsMarkovian
            l{i}=['\lambda = ',num2str(lambda(i))];
        else
            l{i}=['PLR = ',num2str(pgb(i))];
        end
    end
    
    figure(f1); if fourSeq; subplot(2,2,j); hold all; box on; end
    plot(bw/1e3,(NQQ(1:numChains,:).*NQT(1:numChains,:))');
    plot(bw/1e3,max( mnqq( vs.q0 .* ( (bw./vs.R0).*(15/vs.fmax)^-vs.beta_f ).^(-1/vs.beta_q), vs.alpha_q, vs.qmin )...
    .* mnqt( 15, vs.alpha_f, vs.fmax ), mnqq( vs.q0 .* ( (bw./vs.R0) ).^(-1/vs.beta_q), vs.alpha_q, vs.qmin ) ))
    xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); ylim([0.2 1]); 
    if j == 3 || ~fourSeq; legend([l,'PLR = 0'],'Location','Best'); end
    
    if willPlotNQQandNQT
        figure(f2); subplot(2,2,j);
        plot(bw/1e3,NQQ(1:numChains,:)'); ylabel('NQQ');
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(l,'Location',...
                'Best'); end

        figure(f3); subplot(2,2,j);
        plot(bw/1e3,NQT(1:numChains,:)'); ylabel('NQT');
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(l,'Location',...
                'Best'); end 
    end
    
    figure(f4); if fourSeq; subplot(2,2,j); end
    plot(bw/1e3,TotalFECPerc);
    xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); ylim([0 100]); 
    if j == 3 || ~fourSeq; legend(l,'Location','Best'); end

    figure(f5); if fourSeq; subplot(2,2,j); hold all; box on; end
    plot(bw/1e3,F(1:numChains,:)');
    plot(bw/1e3, 15+15*(mnqq( vs.q0 .* ( (bw./vs.R0).*(15/vs.fmax)^-vs.beta_f ).^(-1/vs.beta_q), vs.alpha_q, vs.qmin )...
    .* mnqt( 15, vs.alpha_f, vs.fmax ) < mnqq( vs.q0 .* ( (bw./vs.R0) ).^(-1/vs.beta_q), vs.alpha_q, vs.qmin )) )
    xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); ylim([0 31]); 
    if j == 3 || ~fourSeq; legend([l,'PLR = 0'],'Location','Best'); end
    

    
    meanFECrates = zeros(numChains,numCapacs,L+1);
    for v = 1:numChains
        for u = 1:numCapacs
            vs.f = F(v,u);
            k = ceil( estimFrameSz( vs, R(v,u) ) / PACKET_SIZE );
            m = M{v,u};
            meanFECrates(v,u,:) = 100*find_mean_per_layer( m./(k+m), length(m), L );
        end
    end
    meanFECrates = squeeze( mean(meanFECrates,2) );

    Legend = ['I-Frame'; cell(L,1); 'Overall'];
    for u = 1:L; Legend{u+1} = ['TL(',num2str(u),')']; end    
    
    figure(f6); if fourSeq; subplot(2,2,j); end
    hold all; box on;
    affineModelPlots = zeros(L+2,1);
    if thisIsMarkovian
        horAxis = lambda; horAxisLabel = 'Mean Loss Burst Length (packets)';
        affineModelPlots(L+2) = plot(horAxis,mean(TotalFECPerc));
%         set(gca,'xscale','log');
    else
        horAxis = 100*epsilon; horAxisLabel = 'Packet Loss Rate (%)'; 
        model = fit(horAxis',mean(TotalFECPerc)','poly1');
        fprintf([video,': %fx+%f\n'],model.p1,model.p2);
        x = linspace(min(horAxis),max(horAxis),200);
        affineModelPlots(L+2) = plot(x,x*model.p1+model.p2);
        scatter(horAxis,mean(TotalFECPerc));
    end
    for v = 1:L+1; affineModelPlots(v) = plot(horAxis,meanFECrates(:,v)); end
    xlabel(horAxisLabel); title(video);
    if j == 3 || ~fourSeq; legend( affineModelPlots, Legend, 'Location', 'Best' ); end
    
    if willPlotMonteCarlo && exist([video,'-',num2str(L),'-MC.mat'],'file')
        load([video,'-',num2str(L),'-MC.mat'])

        figure(f7); subplot(2,2,j);
        plot(bw/1e3,30*MeanFrameIntervals(1:numChains,:)'); ylabel('Mean Frame Interval'); ylim([0 15])
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(l,'Location',...
                'Best'); end

        figure(f8); subplot(2,2,j);
        plot(bw/1e3,30*StdFrameIntervals(1:numChains,:)'); ylabel('Std Frame Interval'); ylim([0 10])
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(l,'Location',...
                'Best'); end

        figure(f9); subplot(2,2,j);
        plot(bw/1e3,numFreezes(1:numChains,:)'/100e3); ylabel('Pr(Frz)'); ylim([0 0.025])
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(l,'Location',...
                'Best'); end

        figure(f10); subplot(2,2,j);
        plot(bw/1e3,MeanNumDecFrames(1:numChains,:)'); ylabel('MNDF');
        xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(l,'Location',...
                'Best'); end
    end
    
end

if willSavePlotsAsEps
    target = '~/Google Drive/NYU/Research/papers/fec/fig/';
    saveTightFigure(f1,[target,'quality-',num2str(L),channel,'.eps'])
    saveTightFigure(f4,[target,'videoBitrate-',num2str(L),channel,'.eps'])
    saveTightFigure(f5,[target,'encFrRate-',num2str(L),channel,'.eps'])
    saveTightFigure(f6,[target,'fecRatesPerLayer-',num2str(L),channel,'.eps'])
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