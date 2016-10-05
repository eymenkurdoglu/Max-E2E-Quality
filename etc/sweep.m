close all
file = 'swept.mat';
RBR_IPR = 16/15;
target = '~/Google Drive/NYU/Research/papers/fec/fig/';
if ~exist( file, 'file' )

    VIDEO = 'CREW';
    L = 3;
    PLRs = 0 : 0.05 : 0.2;
    send_rates = 100 : 50 : 1300;
    RESULTS = cell( length(PLRs), length(send_rates) );
    i = 1;
    for PLR = PLRs
        fprintf('# Simulating for PLR = %f\n', PLR);
        j = 1;
        for SR = send_rates
            fprintf('######### Sending rate = %d kbps\n', SR);
            [Qmax, FRopt, RVopt] = heuristic_search( PLR, SR, L, VIDEO );
%             Results = heuristic_search( PLR, SR, L, VIDEO );
            RESULTS{i,j} = Results;
            j = j+1;
        end
        i = i+1;
    end

    save(file, 'RESULTS', 'VIDEO', 'PLRs', 'send_rates', 'L')
else
    load( file )
end

optimalQ = zeros( length(PLRs), length(send_rates) );
optimaleFR = zeros( length(PLRs), length(send_rates) );
optimaldFR = zeros( length(PLRs), length(send_rates) );
optimalRV = zeros( length(PLRs), length(send_rates) );

for i = 1:length(PLRs)
    for j = 1:length(send_rates)
        temp = RESULTS{i,j};
        FR30 = temp{1,1}; FR15 = temp{1,2};
        Q_FR30 = FR30(1,:); Q_FR15 = FR15(1,:);
        if max(Q_FR30) > max(Q_FR15)
            optimalQ(i,j) = max(Q_FR30);
            optimaleFR(i,j) = 30;
            optimaldFR(i,j) = FR30(3,Q_FR30==max(Q_FR30))/RBR_IPR;
            optimalRV(i,j) = FR30(2,Q_FR30==max(Q_FR30));
            if 0
            optimalQ(i,j) = max(Q_FR15);
            optimaleFR(i,j) = 15;
            optimaldFR(i,j) = FR30(3,Q_FR30==max(Q_FR30))/RBR_IPR;
            optimalRV(i,j) = FR15(2,Q_FR15==max(Q_FR15));            
            end
        else
            optimalQ(i,j) = max(Q_FR15);
            optimaleFR(i,j) = 15;
            optimaldFR(i,j) = FR15(3,Q_FR15==max(Q_FR15))/RBR_IPR;
            optimalRV(i,j) = FR15(2,Q_FR15==max(Q_FR15));
        end
    end
end

optimalQ = optimalQ/max(max(optimalQ));
normqoptimalQ = optimalQ ./ repmat(optimalQ(1,:),size(optimalQ,1),1);
figure;
for i = 1:length(PLRs)
    subplot(1,2,1); hold all; box on;
    xlim([min(send_rates),max(send_rates)]/1000); ylabel('Achieved Perceptual Quality'); xlabel('sending rate (Mbps)'); title('IID');
    plot( send_rates/1000, optimalQ(i,:) );
    subplot(1,2,2); hold all; xlim([min(send_rates),max(send_rates)]/1000); ylabel('Achieved Quality Perc. wrt No Loss Case'); xlabel('sending rate (Mbps)'); title('IID');
    plot( send_rates/1000, 100*normqoptimalQ(i,:) );
end
legend('PLR=0','PLR=0.05','PLR=0.10','PLR=0.15','PLR=0.20','Location','SouthEast')
% saveTightFigure([target,'BestQuality-IID.eps'])

% figure; hold all; box on;
% for i = 1:length(PLRs)
%     plot( send_rates/1000, optimaldFR(i,:) );
% end

% figure;
% for i = 1:length(PLRs)
%     subplot(1,2,1); hold all; box on;
%     xlim([min(send_rates),max(send_rates)]/1000); ylim([14 31]);
%     ylabel('Best frame rate (Hz)'); xlabel('sending rate (Mbps)'); title('IID');
%     plot( send_rates/1000, optimaleFR(i,:) );
%     subplot(1,2,2); hold all; box on;
%     xlim([min(send_rates)/1000,max(send_rates)/1000]); ylabel('Best video bitrate (Mbps)'); xlabel('sending rate (Mbps)'); title('IID')
%     plot( send_rates/1000, optimalRV(i,:)/1000 ); 
% end
% legend('PLR=0','PLR=0.05','PLR=0.10','PLR=0.15','PLR=0.20','Location','SouthEast')
% saveTightFigure([target,'BestVRandFR-IID.eps'])

% q30 = RESULTS{5,length(send_rates)}{1,1}(1,:);
% q15 = RESULTS{5,length(send_rates)}{1,2}(1,:);
% video_rates = 50:max(send_rates);
% figure; hold all;
% plot( video_rates(q30>0), q30(q30>0) );
% plot( video_rates(q15>0), q15(q15>0) );