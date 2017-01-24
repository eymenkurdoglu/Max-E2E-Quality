function plotAll( sequences, L )
close all

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;

for j = 1 : length(sequences)
    
    video = sequences{j};
    
    load([video,'-',num2str(L),'.mat'])
    
    numChains = length(pgb);
    numCapacs = length(bw);
    xl = [bw(1) bw(numCapacs)]/1e3;
    
    l = cell(1,numChains); % legend strings
    for i = 1:numChains
        l{i}=['PLR = ',num2str(pgb(i))];
    end
    
    figure(f1); subplot(2,2,j);
    plot(bw/1e3,(NQQ.*NQT)');
    xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); ylim([0 1]); if j == 3; legend(l,'Location','SouthEast'); end
    
    figure(f2); subplot(2,2,j); TotalFECPerc = 100*(1-(R./repmat(bw,numChains,1))');
    plot(bw/1e3,TotalFECPerc);
    xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); ylim([0 100])

    figure(f3); subplot(2,2,j);
    plot(bw/1e3,F');
    xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); ylim([0 31]); if j == 3; legend(l,'Location','SouthEast'); end
    
    Legend = cell(L+2,1);
    Legend{1} = 'I-Frame';
    for u = 1:L
        Legend{u+1} = ['TL(',num2str(u),')'];
    end
    Legend{L+2} = 'Overall';
    
    meanFECrates = zeros(numChains,numCapacs,L+1);
    for v = 1:numChains
        for u = 1:numCapacs
            N = length( M{v,u} );
            vs.f = F(v,u);
            k = ceil( estimFrameSz( vs, R(v,u) ) / PACKET_SIZE );
            m = M{v,u};
            meanFECrates(v,u,:) = 100*find_mean_per_layer( m./(k+m), N, L );
        end
    end
    meanFECrates = squeeze( mean(meanFECrates(:,:,:),2) );
    
    figure(f4); subplot(2,2,j); hold all; box on; title(video); 
    TakeAwayPlots = zeros(L+2,1);
    for v = 1:L+1
        TakeAwayPlots(v) = plot(100*pgb,meanFECrates(:,v));
    end
    model = fit(100*pgb',mean(TotalFECPerc(11:end,:))','poly1');
    fprintf([video,': %fx+%f\n'],model.p1,model.p2);
    scatter(100*pgb,mean(TotalFECPerc(11:end,:)))
    x = linspace(100*min(pgb),100*max(pgb),200); y = x*model.p1+model.p2;   
    TakeAwayPlots(L+2) = plot(x,y);
    xlabel('Packet Loss Rate (%)'); ylabel('FEC bitrate %'); title(video); ylim([0 55])
    if j == 3
        legend(TakeAwayPlots,Legend,'Location','Best');
    end    
    
end

target = '~/Google Drive/NYU/Research/papers/fec/fig/';
saveTightFigure(f1,[target,'quality-',num2str(L),'-IID.eps'])
saveTightFigure(f2,[target,'videoBitrate-',num2str(L),'-IID.eps'])
saveTightFigure(f3,[target,'encFrRate-',num2str(L),'-IID.eps'])
saveTightFigure(f4,[target,'fecRatesPerLayer-',num2str(L),'-IID.eps'])

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