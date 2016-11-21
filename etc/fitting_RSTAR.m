clc
videos = {'CREW','CITY','FOREMAN','HARBOUR','ICE','FG'};
frame_rates = {'30.0','15.0'};
L = 1;
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
    elseif  strcmp(video, 'FG')
        alpha_q =  10.68; alpha_t = 2.8;
    end
 
    q_mean = zeros(2,23);
    MeanBitrate = zeros(2,23);
    
    for fr = 1:length(frame_rates)
        [Lengths, QPs, ~] = parse_log_files( ['~/Dropbox/Matlab/3_OptimalQuality/data/'...
            ,video,'-352x288-',frame_rates{fr},'-',num2str(ipr*str2double(frame_rates{fr}))...
            ,'-',num2str(L),'/'] );
        MeanBitrate(fr,:) = sum(Lengths)*0.008/10;
        q_mean(fr,:) = mean( 2.^((QPs-4)./6) );
    end
    
    Q = (1-exp(-1*alpha_q*(min(min(q_mean))./q_mean)))/(1-exp(-1*alpha_q));
    Q(2,:) = Q(2,:) .* (1-exp(-1*alpha_t*(15./tmax).^0.63))/(1-exp(-1*alpha_t));
    
    subplot(2,length(videos),v); hold all; title(video);
    plot( MeanBitrate(1,:), Q(1,:) ); xlim([0 1200]); ylim([0.2 1])
    plot( MeanBitrate(2,:), Q(2,:) )
    xlabel('bitrate')
    ylabel('NQQ*NQT')
    legend(frame_rates,'Location','SouthEast')
    
    R_max = max(MeanBitrate(1,:));
    q_min = min(q_mean(1,:));
    
%     [model, ~, ~] = fit((log(q_min./q_mean(1,:)))', (log(MeanBitrate(1,:)/R_max))', 'poly1');  
%     fprintf('a_q=%f\n',model.p1);
    
    [model, ~, ~] = fit((log(q_min./q_mean(2,:)))', (log(MeanBitrate(2,:)/R_max))', 'poly1');  
    fprintf([videos{v},': beta_q = %f; beta_t = %f; Rmax = %f; qmin = %f;\n'],model.p1,-1*(model.p2)/log(2),R_max,q_min);
    
    subplot(2,length(videos),v+length(videos)); hold all;
    loglog( q_min./q_mean(2,:), MeanBitrate(2,:)/R_max );
    xlabel('Normalized QS')
    ylabel('Normalized bitrate')
end