function frameSizeModelPlot( sequences, L )

dbstop if error
close all

videoStates = cell(size(sequences));
intraperiodDur = 16/15;

handle = figure; hold all; box on;
xlabel('Distance to reference frame (msec)'); ylabel('Normalized mean P-frame size');

lineArray = [];%cell(2*length(sequences),1);
legenD = cell(2*length(sequences),1);

i = 1;
for video = sequences
    
     vs = struct( 'video', video{1}, 'ipr', intraperiodDur, 'vidsize', '704x576', ...
        'fr', [30 15], 'L', L, 'maxValidTargetBitrate', 1300 );    
    
    vs.eta = zeros( size(vs.fr) );
    vs.theta = zeros( size(vs.fr) );
    vs.fitPlots = cell( fliplr(size(vs.fr)) );
    
    videoStates{i} = frameSizeModel( vs, handle );
    legenD{2*i-1} = [vs.video,'-15 Hz']; legenD{2*i} = [vs.video,'-30 Hz'];
    lineArray = [lineArray; videoStates{i}.fitPlots{2}];
    lineArray = [lineArray; videoStates{i}.fitPlots{1}];
%     lineArray{2*i-1} = videoStates{i}.fitPlots{2}; lineArray{2*i} = videoStates{i}.fitPlots{1};
    i = i+1;
end

legend( lineArray, legenD, 'Location', 'SouthEast' )
saveTightFigure(handle,'~/Google Drive/NYU/Research/papers/fec/fig/frszmodel.eps')
return

function vs = frameSizeModel( vs, varargin )

    for f = vs.fr
    
        path = ['data/',vs.video,'-',vs.vidsize,'-',strcat(...
            num2str(f),'.0'),'-',num2str(f * vs.ipr),'-',num2str(vs.L),'/'];
        [Lengths, ~, Target_BitRates] = parse_log_files( path ); 

        % choose valid bitrates if need be
        valid = Target_BitRates <= vs.maxValidTargetBitrate;
        numValid = sum(valid);
        Lengths = Lengths(:,valid);

        meanFrmSz = zeros( vs.L+1, numValid ); % average frame lengths in each layer per video bitrate

        for r = 1:numValid
            meanFrmSz(:,r) = find_mean_per_layer( Lengths(:,r), f * vs.ipr, vs.L ); % BE CAREFUL!!!!!
        end

        meanFrmSz = meanFrmSz./repmat(meanFrmSz(1,:),vs.L+1,1); % mean P-frame length/mean I-frame length
        meanFrmSz(1,:) = [];

        model = fit([1000*fliplr(2.^(0:vs.L-1))/f, 0]', [mean(meanFrmSz,2); 0],...
            fittype('fitmodel( x, a, b )'), 'Startpoint', [0 0]);
        
        if nargin == 2
            x = 1000*fliplr(2.^(0:vs.L-1))/f;
            y = mean(meanFrmSz,2);
            figure(varargin{1});
            if f == 30; line = '-'; else line = '--'; end            
            if strcmp(vs.video,'CREW'); color = 'r';
            elseif strcmp(vs.video,'CITY'); color = 'g';
            elseif strcmp(vs.video,'HARBOUR'); color = 'b';
            elseif strcmp(vs.video,'SOCCER'); color = 'c';
            end
            scatter( x, y, color )
            vs.fitPlots{ f == vs.fr } = plot(linspace(0,max(x),200),fitmodel( linspace(0,max(x),200)...
                , model.a, model.b ),[color,line],'LineWidth',2);
            
        end
        
        vs.theta( f == vs.fr ) = model.a;
        vs.eta  ( f == vs.fr ) = model.b;  
        
        fprintf('%s@%d Hz => theta=%f, eta=%f\n',vs.video,  f, model.a,  model.b );
    
    end
    
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