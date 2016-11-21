close all
clc

target = '~/Google Drive/NYU/Research/papers/fec/fig/';
RBR_FRAME_RATES = {'30.0','15.0'};
RBR_NUM_TEMP_LYRS = 1;
RBR_VIDEOS = {'CREW','CITY','FOREMAN','HARBOUR','ICE','FG'};
RBR_IPR = 16/15;

frszmodel = figure; hold all; box on;
xlabel('Distance to reference frame (msec)'); ylabel('Normalized mean P-frame size');

% frszwrtrate = figure;
% xlabel('Target bitrate (kbps)'); ylabel('Normalized P-frame size'); 

leg = {};
YO = [];

for fr = 1:length(RBR_FRAME_RATES)

f = str2double(RBR_FRAME_RATES{fr});

N = RBR_IPR * f;

if f == 30
    line = '-';
else
    line = '--';
end

    for v = 1:length(RBR_VIDEOS)

        if strcmp(RBR_VIDEOS{v},'CREW')
            color = 'r';
        elseif strcmp(RBR_VIDEOS{v},'CITY')
            color = 'm';
        elseif strcmp(RBR_VIDEOS{v},'FOREMAN')
            color = 'b';
        elseif strcmp(RBR_VIDEOS{v},'HARBOUR')
            color = 'k';
        elseif strcmp(RBR_VIDEOS{v},'ICE')
            color = 'g';
        elseif strcmp(RBR_VIDEOS{v},'FG')
            color = 'c';
        end

        [Lengths, QPs, Target_BitRates] = parse_log_files(['~/Dropbox/Matlab/3_OptimalQuality/data/'...
            ,RBR_VIDEOS{v},'-352x288-',RBR_FRAME_RATES{fr},'-',num2str(N),'-',num2str(RBR_NUM_TEMP_LYRS),'/']); 

        avgfrmlen = zeros(RBR_NUM_TEMP_LYRS+1,length(Target_BitRates)); % average frame lengths in each layer per video bitrate

        for r = 1:length(Target_BitRates)
            avgfrmlen(:,r) = find_mean_per_layer( Lengths(:,r), N, RBR_NUM_TEMP_LYRS );
        end

        avgfrmlen = avgfrmlen./repmat(avgfrmlen(1,:),RBR_NUM_TEMP_LYRS+1,1); % P-frame length/I-frame length
        avgfrmlen(1,:) = [];
        
        mean(avgfrmlen,2)
        std(avgfrmlen')./mean(avgfrmlen')

%         figure(frszwrtrate); 
%         subplot(length(RBR_VIDEOS), length(RBR_FRAME_RATES), (v-1)*length(RBR_FRAME_RATES)+fr, 'Parent',frszwrtrate,'XTick',[0.2 1.2]);
%         hold all; box on; title(RBR_VIDEOS{v});
%         xlim([0.150, 1.200]); ylim([0 1])
%         for l = 1:RBR_NUM_TEMP_LYRS
%             plot( Target_BitRates/1000, avgfrmlen(l,:) )
%         end
%         if strcmp( RBR_VIDEOS{v}, 'FG' )
%             xlabel('Video bitrate (Mbps)')
%         end

        x = [1000*fliplr(2.^(0:RBR_NUM_TEMP_LYRS-1))/f, 0]';
        y = [mean(avgfrmlen,2); 0];

        ft = fittype('fitmodel( x, a, b )');
        model = fit(x, y, ft, 'Startpoint', [0 0]);    
        fprintf([RBR_VIDEOS{v},', FR=%d, a=%f, b=%f\n'],f,model.a,model.b);

        figure(frszmodel);
        scatter( x, y, color );
        HOBELE = plot(linspace(0,max(x),200),fitmodel( linspace(0,max(x),200), model.a, model.b ),[color,line],'LineWidth',2);
        YO = [YO, HOBELE];
        leg = [leg; '', [RBR_VIDEOS{v},'-',RBR_FRAME_RATES{fr},' Hz']];
    end
end
legend(YO,leg,'Location','SouthEast')
% saveTightFigure(frszwrtrate,[target,'frszwrtrate.eps'])
% saveTightFigure(frszmodel,[target,'frszmodel.eps'])