dbstop if error
clear all
clc

global alpha_q minQS alpha_t maxt ipr

video = 'CREW';
alpha_q = 4.51; minQS = 1000; alpha_t = 3.09; maxt = 30;
ipr = 16/15;
frame_rates = {'30.0','15.0'};
templayers = 3;

RBR_PACKET_SIZE = 200;
RBR_SELECTED_IPR = 4;

Packets = cell(1,length(frame_rates));
QS = cell(1,length(frame_rates));

for fr = 1:length(frame_rates)
    N = ipr*str2double(frame_rates{fr});
    path = ['~/Dropbox/Matlab/3_OptimalQuality/data/',video,'-352x288-',frame_rates{fr},'-',num2str(N),'-3/'];

    Lengths = []; QPs = []; target_bitrates = [];
    for file = dir( strcat(path,'*.txt') )'
        contents = dlmread( strcat(path,file.name) );
        Lengths  = [Lengths contents(:,1)];
        QPs = [QPs contents(:,2)];
        dum = strsplit(file.name,'.');
        target_bitrates = [target_bitrates str2double(dum{1})];
    end
    
    % sort in ascending target bitrates instead of alphabetical
    [target_bitrates,ix] = sort(target_bitrates); Lengths = Lengths(:,ix); QPs = QPs(:,ix);
    
    % convert QPs to QSs
    QSs = 2.^((QPs-4)./6);
    
    % packetize the frames in the selected intra-period
    range = ((RBR_SELECTED_IPR-1)*N + 1) : RBR_SELECTED_IPR*N;
    Packets{fr} = ceil( Lengths(range,:) / RBR_PACKET_SIZE );
    
    QS{fr} = mean(QSs(range,:));
    minQS = min( minQS, min(QS{fr}) );
end  

alphas = .05 : .05 : .4;
betas = 1-alphas;
abws = round( (100 : 50 : 1200)/8e-3/200 );

results = cell(length(alphas),length(abws));
for i = 1:length(alphas)
    alpha = alphas(i);
%     for beta = betas
        for j = 1:length(abws)
            abw = abws(j);
            res = zeros(length(frame_rates),4);
%             best_qual = 0; best_fr = 0; best_sr = 0; best_num_frames = 0;
            for fr = 1:length(frame_rates)
                fec = abw - sum( Packets{fr} ); % some will be negative
                qp = 6*log2(QS{fr})+4;
                [m, quality, num_frames] = find_max_qual( alpha, 1-alpha, fec, Packets{fr}, QS{fr}, templayers );
                if ~isempty(m)
                    res(fr,:) = [quality, 100*(abw-fec(m))/abw, num_frames, qp(m)];
                end
%                 if quality > best_qual
%                     best_fr = fr; best_sr = abw-fec(m); best_qual = quality; best_num_frames = num_frames;
%                 end
            end
            results{i,j} = res;
%             results{i,j} = [best_fr, best_sr, best_qual, best_num_frames];
        end
%     end
end