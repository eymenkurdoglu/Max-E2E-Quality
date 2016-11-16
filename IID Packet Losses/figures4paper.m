RBR_IPR = 16/15;
RBR_PACKET_SIZE = 200;

close all
target = '~/Google Drive/NYU/Research/papers/fec/fig/';
f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;

videos = {'CREW','CITY','HARBOUR','FG'};

for j = 1:length(videos)
    load(['optimized',videos{j},'.mat'])
    l=cell(1,length(PLRs));
    for i = 1:length(PLRs)
        l{i}=['PLR=',num2str(PLRs(i))];
    end
    SRs = SRs/1000;
    optimalRV = optimalRV/1000;
    %###############################################################
    figure(f1); subplot(2,2,j);
    plot(SRs,optimalQ'/max(max(optimalQ))); xlim([SRs(1) SRs(end)]);
    % xlabel('Sending Bitrate (Mbps)')
    title(videos{j}); 
    if j == 4 
        legend(l,'Location','SouthEast') 
    end
    %###############################################################
    figure(f2); subplot(2,2,j);
    plot(SRs,100*(1-(optimalRV./repmat(SRs,length(PLRs),1))')); xlim([SRs(1) SRs(end)]); ylim([10 100]);
    % xlabel('Sending Bitrate (Mbps)'); ylabel('%')
    title(videos{j}); 
    if j == 1
        legend(l,'Location','Best')
    end
    %###############################################################
    figure(f3); subplot(2,2,j);
    plot(SRs,optimaleFR'); xlim([SRs(1) SRs(end)]); ylim([0 31]);
    % xlabel('Sending Bitrate (Mbps)'); ylabel('Encoding Frame Rate (Hz)')
    title(videos{j}); 
    if j == 4
        legend(l,'Location','SouthEast');
    end
    %###############################################################
    optimaldFR = zeros(length(PLRs),length(SRs));
    for v = 1:length(PLRs)
        for u = 1:length(SRs)
            optimaldFR(v,u) = optimalPMF{v,u} * (0:length(optimalPMF{v,u})-1)';
        end
    end
    optimaldFR = optimaldFR / (16/15);

    figure(f4); subplot(2,2,j);
    plot(SRs,optimaldFR'); xlim([SRs(1) SRs(end)]); ylim([0 31]);
    % xlabel('Sending Bitrate (Mbps)'); ylabel('Hz');
    title(videos{j});
    if j == 4
        legend(l,'Location','SouthEast');
    end
    %###############################################################
    FECrate = zeros(length(PLRs),length(SRs),L+1);
    for v = 1:length(PLRs)
        for u = 1:length(SRs)
            N = length(optimalm{v,u});
            k = ceil( estimate_frame_lengths(optimaleFR(v,u), L, N, RBR_IPR*1e3*optimalRV(v,u)/0.008, videos{j}) / RBR_PACKET_SIZE );
            m = optimalm{v,u};
            r = find_mean_per_layer( k./(k+m), N, L );
            FECrate(v,u,:) = r;
        end
    end
    figure(f5); subplot(2,2,j);
    plot(PLRs,squeeze(mean(FECrate(:,:,:),2))); xlim([PLRs(1) PLRs(end)]); %ylim([0 0.6])
    % xlabel('Packet Loss Rate');
    title(videos{j}); 
    if j == 2
        legend('I-frame','TL(1)','TL(2)','TL(3)','Location','Best');
    end
end
saveTightFigure(f1,[target,'BestQuality-IID.eps'])
saveTightFigure(f2,[target,'BestVideoBitrate-IID.eps'])
saveTightFigure(f3,[target,'BestEncFrRate-IID.eps'])
saveTightFigure(f4,[target,'BestDecFrRate-IID.eps'])
saveTightFigure(f5,[target,'FECRatesPerLayer-IID.eps'])