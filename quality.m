% close all
clear all
clc
videos = {'CREW','CITY','FOREMAN','HARBOUR'};
frame_rates = {'30.0','15.0'};
ipr = 16/15;
figure;
for v = 1:length(videos)

    video = videos{v};
    
    tmax = 30;

    if      strcmp(video, 'CREW')
        alpha_q = 4.51; alpha_t = 3.09; 
    elseif  strcmp(video, 'CITY')
        alpha_q = 7.25; alpha_t = 4.10; 
    elseif  strcmp(video, 'HARBOUR')
        alpha_q = 9.65; alpha_t = 2.83; 
    elseif  strcmp(video, 'ICE')    
        alpha_q = 5.61; alpha_t = 3.00; 
    elseif  strcmp(video, 'FOREMAN')    
        alpha_q = 4.57; alpha_t = 3.80;
    end
    
    subplot(2,4,v); hold all;
    title(video)
    q_mean = zeros(2,23);
    qp_mean = zeros(2,23);
    bitrate_mean = zeros(2,23);
    for fr = 1:length(frame_rates)
        N = ipr*str2double(frame_rates{fr});
        path = ['~/Dropbox/Matlab/3_OptimalQuality/data/',video,'-352x288-',frame_rates{fr},'-',num2str(N),'-3/'];

        lengths = []; qp = []; rates = [];
        for file = dir( strcat(path,'*.txt') )'
            contents = dlmread( strcat(path,file.name) );
            lengths  = [lengths contents(:,1)];
            qp = [qp contents(:,2)];
            dum = strsplit(file.name,'.');
            rates = [rates str2double(dum{1})];
        end
        [rates,ix] = sort(rates);
        lengths = lengths(:,ix);
%         sum(lengths(:,ix))*0.008/10
        qp = qp(:,ix);
        q = 2.^((qp-4)./6);
        q_mean(fr,:) = mean(q);
        qp_mean(fr,:) = mean(qp);
        bitrate_mean(fr,:) = sum(lengths)*0.008/10;
    end  
    
%     q = 2.^((qp_mean-4)./6); % QP = 4+6*log2(q)
    qmin = min(min(q_mean));
    Q = (1-exp(-1*alpha_q*(qmin./q_mean)))/(1-exp(-1*alpha_q));
    Q(2,:) = Q(2,:) .* (1-exp(-1*alpha_t*(15./tmax).^0.63))/(1-exp(-1*alpha_t));
    plot( bitrate_mean(1,:), Q(1,:) ); xlim([0 1200]); ylim([0.2 1])
    plot( bitrate_mean(2,:), Q(2,:) )
    xlabel('bitrate')
    ylabel('NQQ*NQT')
    legend(frame_rates,'Location','SouthEast')
    
    subplot(2,4,v+4); hold all;
    plot( bitrate_mean(1,:), qp_mean(1,:) ); xlim([0 1200])
    plot( bitrate_mean(2,:), qp_mean(2,:) ); ylim([20 50])
    ylabel('QP')
    xlabel('bitrate')
end