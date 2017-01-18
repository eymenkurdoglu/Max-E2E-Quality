function plotAll( sequences, L )
close all
% target = '~/Google Drive/NYU/Research/papers/fec/fig/';

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;

for j = 1 : length(sequences)
    
    video = sequences{j};
    
    load([video,'-',num2str(L),'.mat'])
    
    numChains = length(pgb);
    numCapacs = length(bw);
    
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
    
%     bw = bw/1000; R = R/1000;
    
    figure(f1); subplot(2,2,j);
    plot(bw,(NQQ.*NQT)');
    xlabel('Sending Bitrate (Mbps)'); title(video); 
    
    figure(f2); subplot(2,2,j); TotalFECPerc = 100*(1-(R./repmat(bw,numChains,1))');
    plot(bw,TotalFECPerc);
    xlabel('Sending Bitrate (Mbps)'); ylabel('%')
    title(video); 

    figure(f3); subplot(2,2,j);
    plot(bw,F'); ylim([0 31]);
    xlabel('Sending Bitrate (Mbps)'); ylabel('Encoding Frame Rate (Hz)')
    title(video);
    
    Legend = cell(L+2,1);
    Legend{1} = 'I-Frame';
    for u = 1:L
        Legend{u+1} = ['TL(',num2str(u),')'];
    end
    Legend{L+2} = 'Overall';
    
    if ~isfield(vs,'eta')
        eta = zeros( size(vs.fr) );
        theta = eta;
        for f = vs.fr
            vs.f = f;
            vs = frameSizeModel( vs );
            eta( vs.fr == f ) = vs.eta;
            theta( vs.fr == f ) = vs.theta;
        end
        vs.eta = eta; vs.theta = theta;
        save([video,'-',num2str(L),'.mat'],'NQQ','NQT','F','R','M','D','pgb','pbg','bw','numLayers','vs','PACKET_SIZE')
    end
    
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
    meanFECrates = squeeze( mean(meanFECrates(:,11:end,:),2) );
    
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
    if j == 3
        legend(TakeAwayPlots,Legend,'Location','Best');
    end    
    
end

% saveTightFigure(f1,[target,'BestQuality-IID.eps'])
% saveTightFigure(f2,[target,'BestVideoBitrate-IID.eps'])
% saveTightFigure(f3,[target,'BestEncFrRate-IID.eps'])
% saveTightFigure(f4,[target,'BestDecFrRate-IID-IPP.eps'])
% saveTightFigure(f5,[target,'FECRatesPerLayer-IID.eps'])
% saveTightFigure(f6,[target,'FECRatesTakeAway-IID-IPP.eps'])

return

function vs = frameSizeModel( vs )
    
    maxValidTargetBitrate = 1200;
    
    path = ['data/',vs.video,'-',vs.vidsize,'-',strcat(...
        num2str(vs.f),'.0'),'-',num2str(vs.f * vs.ipr),'-',num2str(vs.L),'/'];
    [Lengths, ~, Target_BitRates] = parse_log_files( path ); 
    
    % choose valid bitrates if need be
    valid = Target_BitRates <= maxValidTargetBitrate;
    numValid = sum(valid);
    Lengths = Lengths(:,valid);
    
    meanFrmSz = zeros( vs.L+1, numValid ); % average frame lengths in each layer per video bitrate

    for r = 1:numValid
        meanFrmSz(:,r) = find_mean_per_layer( Lengths(:,r), vs.f * vs.ipr, vs.L ); % BE CAREFUL!!!!!
    end

    meanFrmSz = meanFrmSz./repmat(meanFrmSz(1,:),vs.L+1,1); % mean P-frame length/mean I-frame length
    meanFrmSz(1,:) = [];

    model = fit([1000*fliplr(2.^(0:vs.L-1))/vs.f, 0]', [mean(meanFrmSz,2); 0],...
        fittype('fitmodel( x, a, b )'), 'Startpoint', [0 0]);

    vs.theta = model.a;
    vs.eta = model.b;
    
    fprintf('%s@%d Hz => theta=%f, eta=%f\n',vs.video,  vs.f, vs.theta,  vs.eta);
    
return

function k = estimFrameSz( vs, R )
% estimFrameSz      estimate the sizes of the I and P frames in the
% intra-period given the target encoding bitrate
%  Size of the P-frame tau ms apart from its reference normalized wrt the
%  I-frame is estimated using the model: | zhat | = 1 - exp(-theta*tau^-eta)
%  

    f = vs.f;
    L = vs.L;
    N = vs.f * vs.ipr; % number of frames in intra-period
    B = 1000 * vs.ipr * R/8; % byte budget in intra-period
    
    eta = vs.eta( vs.fr == vs.f );
    theta = vs.theta( vs.fr == vs.f );
    
    zhatNorm = fliplr( 1-exp(-theta*(1000*(2.^(0:L-1))./f).^eta) );

    % n: number of P-frames in intra-period per TL, e.g. [7 8 16]
    n = fliplr( N ./ (2.^[ 1:L-1, L-1]) );
    n(1) = n(1)-1;

    zhat = (B(:)/(1 + n*zhatNorm')*[1, zhatNorm])';

    % find frame indices for each TL and create the k vector or matrix
    % start from the last layer, then previous, ...
    k = zeros(N,length(B));
    for r = 1:length(B)
        for layer = L : -1 : 1   
            k(mod(0:N-1,2^(L-layer))==0,r) = zhat(layer+1,r);
        end
    end
    k(1,:) = zhat(1,:);
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