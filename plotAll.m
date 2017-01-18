function plotAll( sequences )

% target = '~/Google Drive/NYU/Research/papers/fec/fig/';

f1 = figure;
f2 = figure;
f3 = figure;

for j = 1:length(sequences)
    
    video = sequences{j};
    
    load([video,'-3.mat'])
    
    numChains = length(pgb);
    
%     if strcmp(videos{j},'CREW')
%         color = 'r';
%     elseif strcmp(videos{j},'CITY')
%         color = 'm';
%     elseif strcmp(videos{j},'FOREMAN')
%         color = 'b';
%     elseif strcmp(videos{j},'HARBOUR')
%         color = 'k';
%     elseif strcmp(videos{j},'ICE')
%         color = 'g';
%     elseif strcmp(videos{j},'SOCCER')
%         color = 'c';
%     end
    
    l = cell(1,numChains); % legend strings
    for i = 1:numChains
        l{i}=['PLR = ',num2str(pgb(i))];
    end
    
    bw = bw/1000; R = R/1000;
    
    figure(f1); subplot(2,2,j);
    plot(bw,(NQQ.*NQT)');
    xlabel('Sending Bitrate (Mbps)'); title(video); 
    
    figure(f2); subplot(2,2,j);
    plot(bw,100*(1-(R./repmat(bw,numChains,1))'));
    xlabel('Sending Bitrate (Mbps)'); ylabel('%')
    title(video); 

    figure(f3); subplot(2,2,j);
    plot(bw,F'); ylim([0 31]);
    xlabel('Sending Bitrate (Mbps)'); ylabel('Encoding Frame Rate (Hz)')
    title(video);
    
%     Legend = cell(L+2,1);
%     Legend{1} = 'I-Frame';
%     for u = 1:L
%         Legend{u+1} = ['TL(',num2str(u),')'];
%     end
%     Legend{L+2} = 'Overall';
%     
%     meanFECrates = zeros(length(pgb),length(bw),L+1);
%     for v = 1:length(pgb)
%         for u = 1:length(bw)
%             N = length(M{v,u});
%             k = ceil( estimate_frame_lengths(optimaleFR(v,u), L, N, RBR_IPR*1e3*optimalRV(v,u)/0.008, videos{j}) / RBR_PACKET_SIZE );
%             m = M{v,u};
%             meanFECrates(v,u,:) = 100*find_mean_per_layer( m./(k+m), N, L );
%         end
%     end
%     meanFECrates = squeeze( mean(meanFECrates(:,11:end,:),2) );
%     
%     figure(f5); subplot(2,2,j); hold on; box on; title(videos{j}); 
%     xlim([5 25]); ylim([0 60]);
%     TakeAwayPlots = zeros(L+2,1);
%     for v = 1:L+1
%         TakeAwayPlots(v) = plot(100*pgb,meanFECrates(:,v));
%     end
%     model = fit(100*pgb',mean(TotalFECPerc(11:end,:))','poly1');
%     fprintf([videos{j},': %fx+%f\n'],model.p1,model.p2);
%     scatter(100*pgb,mean(TotalFECPerc(11:end,:)))
%     x = linspace(100*min(pgb),100*max(pgb),200); y = x*model.p1+model.p2;   
%     TakeAwayPlots(L+2) = plot(x,y);
%     if j == 3
%         legend(TakeAwayPlots,Legend,'Location','Best');
%     end    
    
end

% saveTightFigure(f1,[target,'BestQuality-IID.eps'])
% saveTightFigure(f2,[target,'BestVideoBitrate-IID.eps'])
% saveTightFigure(f3,[target,'BestEncFrRate-IID.eps'])
% saveTightFigure(f4,[target,'BestDecFrRate-IID-IPP.eps'])
% saveTightFigure(f5,[target,'FECRatesPerLayer-IID.eps'])
% saveTightFigure(f6,[target,'FECRatesTakeAway-IID-IPP.eps'])

return