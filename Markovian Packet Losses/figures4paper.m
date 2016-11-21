RBR_IPR = 16/15;
RBR_PACKET_SIZE = 200;

close all
target = '~/Google Drive/NYU/Research/papers/fec/fig/';
f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;

% videos = {'CREW','CITY','HARBOUR','FG'};
videos = {'CREW'};

for j = 1:length(videos)
    
    load(['optimized',videos{j},'.mat'])
    
    EBL = 1./P_BGs;
    lgnd=cell(1,length(P_GBs));
    for i = 1:length(P_GBs)
        lgnd{i}=['\lambda=',num2str(EBL(i))];
    end
    SRs = SRs/1000;
    optimalRV = optimalRV/1000;
    %###############################################################
    figure(f1); subplot(2,2,j);
    plot(SRs,optimalQ'/max(max(optimalQ))); xlim([SRs(1) SRs(end)]);
    % xlabel('Sending Bitrate (Mbps)')
    title(videos{j}); 
    if j == 4 
        legend(lgnd,'Location','SouthEast') 
    end
    %###############################################################
    figure(f2); subplot(2,2,j);
    plot(SRs,100*(1-(optimalRV./repmat(SRs,length(P_GBs),1))')); xlim([SRs(1) SRs(end)]); ylim([0 100]);
    % xlabel('Sending Bitrate (Mbps)'); ylabel('%')
    title(videos{j}); 
    if j == 1
        legend(lgnd,'Location','Best')
    end
    %###############################################################
    figure(f3); subplot(2,2,j);
    plot(SRs,optimaleFR'); xlim([SRs(1) SRs(end)]); ylim([0 31]);
    % xlabel('Sending Bitrate (Mbps)'); ylabel('Encoding Frame Rate (Hz)')
    title(videos{j}); 
    if j == 4
        legend(lgnd,'Location','SouthEast');
    end
    %###############################################################
    optimaldFR = zeros(length(P_GBs),length(SRs));
    for v = 1:length(P_GBs)
        for u = 1:length(SRs)
            optimaldFR(v,u) = (0:length(optimalPMF{v,u})-1) * optimalPMF{v,u};
        end
    end
    optimaldFR = optimaldFR / (16/15);

    figure(f4); subplot(2,2,j);
    plot(SRs,optimaldFR'); xlim([SRs(1) SRs(end)]); ylim([0 31]);
    % xlabel('Sending Bitrate (Mbps)'); ylabel('Hz');
    title(videos{j});
    if j == 4
        legend(lgnd,'Location','SouthEast');
    end
    %###############################################################
    FECrate = zeros(length(P_GBs),length(SRs),L+1);
    for v = 1:length(P_GBs)
        for u = 1:length(SRs)
            N = length(optimalm{v,u});
            k = ceil( estimate_frame_lengths(optimaleFR(v,u), L, N, RBR_IPR*1e3*optimalRV(v,u)/0.008, videos{j}) / RBR_PACKET_SIZE );
            m = optimalm{v,u};
            r = find_mean_per_layer( k./(k+m), N, L );
            FECrate(v,u,:) = r;
        end
    end
    figure(f5); subplot(2,2,j);
    plot(EBL,squeeze(mean(FECrate(:,:,:),2))); %xlim([P_GBs(1) P_GBs(end)]); %ylim([0 0.6])
    % xlabel('Packet Loss Rate');
    title(videos{j}); 
    if j == 2
        legend('I-frame','TL(1)','TL(2)','TL(3)','Location','Best');
    end
end
% saveTightFigure(f1,[target,'BestQuality-Mark.eps'])
% saveTightFigure(f2,[target,'BestVideoBitrate-Mark.eps'])
% saveTightFigure(f3,[target,'BestEncFrRate-Mark.eps'])
% saveTightFigure(f4,[target,'BestDecFrRate-Mark.eps'])
% saveTightFigure(f5,[target,'FECRatesPerLayer-Mark.eps'])