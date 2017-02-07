function plotDiff( sequences )
% close all

%% 2x2 plots, section a

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;
f6 = figure;
f7 = figure;
f8 = figure;
f9 = figure;

for j = 1 : length(sequences)
    
    video = sequences{j};
    
    for L = [1 3]
    
        load([video,'-',num2str(L),'.mat'])
    
        numChains = length(pgb);
        numCapacs = length(bw);
        xl = [bw(1) bw(numCapacs)]/1e3;
    
        l = cell(1,numChains); % legend strings
        for i = 1:numChains
            l{i}=['PLR = ',num2str(pgb(i))];
        end

        if L == 1; Q = (NQQ.*NQT)'; else Q = Q - (NQQ.*NQT)'; end
        if L == 1; nqq = NQQ'; else nqq = nqq - NQQ'; end
        if L == 1; nqt = NQT'; else nqt = nqt - NQT'; end
        if L == 1; TotalFECPerc = (100*(1-(R./repmat(bw,numChains,1))'))'; else TotalFECPerc =...
                TotalFECPerc - (100*(1-(R./repmat(bw,numChains,1))'))'; end
        load([video,'-',num2str(L),'-MC.mat'])
        if L == 1; MFI = MeanFrameIntervals'; else MFI = MFI - MeanFrameIntervals'; end
        if L == 1; SFI = StdFrameIntervals'; else SFI = SFI - StdFrameIntervals'; end
        if L == 1; FRP = numFreezes'/100e3; else FRP = FRP - numFreezes'/100e3; end
        if L == 1; MNDF = MeanNumDecFrames'; else MNDF = MNDF - MeanNumDecFrames'; end
        if L == 1
            x = 100*(1-(R./repmat(bw,numChains,1))');
            model = fit(100*pgb',mean(x(11:end,:))','poly1');
            y = (linspace(100*min(pgb),100*max(pgb),200))*model.p1+model.p2;
        else
            x = 100*(1-(R./repmat(bw,numChains,1))');
            model = fit(100*pgb',mean(x(11:end,:))','poly1');
            y = y - ((linspace(100*min(pgb),100*max(pgb),200))*model.p1+model.p2);
        end
        fprintf([video,'-',num2str(L),': %fx+%f\n'],model.p1,model.p2);
    end
        
    figure(f1); subplot(2,2,j);
    plot(bw/1e3,Q); if  j == 1; ylabel('Q_{IPP}-Q_{hPP}'); end
    xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(l,'Location',...
            'NorthEast'); end
    
    figure(f4); subplot(2,2,j);
    plot(bw/1e3,TotalFECPerc); 
    if  j == 1; ylabel('FecPerc_{IPP}-FecPerc_{hPP}'); end
    if j == 3; legend(l,'Location','NorthEast'); end
    xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); ylim([0 100]);
    
    figure(f5); subplot(2,2,j);
    plot(bw/1e3,30*MFI); if  j == 1; ylabel('AvgFrInt_{IPP}-AvgFrInt_{hPP}'); end
    xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(l,'Location',...
            'NorthEast'); end
    
    figure(f6); subplot(2,2,j);
    plot(bw/1e3,30*SFI); if  j == 1; ylabel('StdFrInt_{IPP}-StdFrInt_{hPP}'); end
    xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(l,'Location',...
            'NorthEast'); end
    
    figure(f7); subplot(2,2,j);
    plot(bw/1e3,FRP); if  j == 1; ylabel('Pr_{IPP}(Frz)-Pr_{hPP}(Frz)'); end
    xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(l,'Location',...
            'NorthEast'); end
    
    figure(f2); subplot(2,2,j);
    plot(bw/1e3,nqq); if  j == 1; ylabel('NQQ_{IPP}-NQQ_{hPP}'); end
    xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(l,'Location',...
            'NorthEast'); end
    
    figure(f3); subplot(2,2,j);
    plot(bw/1e3,nqt); ylabel('NQT_{IPP}-NQT_{hPP}');
    xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(l,'Location',...
            'SouthEast'); end
    
    figure(f8); subplot(2,2,j);
    plot(bw/1e3,MNDF); ylabel('AvgNumDecFr_{IPP}-AvgNumDecFr_{hPP}');
    xlabel('Sending Bitrate (Mbps)'); title(video); xlim(xl); if j == 3; legend(l,'Location',...
            'SouthEast'); end       
    
    figure(f9); subplot(2,2,j); hold all; box on; title(video);
    plot(linspace(100*min(pgb),100*max(pgb),200),y);
    xlabel('Packet Loss Rate (%)'); ylabel('FecPerc_{IPP}-FecPerc_{hPP}'); title(video);
    
end

% target = '~/Google Drive/NYU/Research/papers/fec/fig/';
% saveTightFigure(f1,[target,'quality-',num2str(L),'-IID.eps'])
% saveTightFigure(f2,[target,'videoBitrate-',num2str(L),'-IID.eps'])
% saveTightFigure(f3,[target,'encFrRate-',num2str(L),'-IID.eps'])
% saveTightFigure(f4,[target,'fecRatesPerLayer-',num2str(L),'-IID.eps'])

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