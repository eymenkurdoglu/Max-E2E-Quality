function vs = initState(video,numLayers,Rmax,intraperiodDur)
% This function initializes the video state object, which holds the Q-STAR
% and R-STAR parameters of the selected video sequence, as well as the
% parameters that summarize the coding structure used.

    % Get Q-STAR parameters
    if      strcmp(video, 'CREW')
        alpha_q = 4.51; alpha_f = 3.09; 
    elseif  strcmp(video, 'CITY')
        alpha_q = 7.25; alpha_f = 4.10; 
    elseif  strcmp(video, 'HARBOUR')
        alpha_q = 9.65; alpha_f = 2.83; 
    elseif  strcmp(video, 'ICE')    
        alpha_q = 5.61; alpha_f = 3.00; 
    elseif  strcmp(video, 'FOREMAN')    
        alpha_q = 4.57; alpha_f = 3.80;
    elseif  strcmp(video, 'FG')
        alpha_q = 10.68; alpha_f = 2.8;
    elseif  strcmp(video, 'SOCCER')
        alpha_q = 6.31; alpha_f = 2.23;
    end

    vs = struct( 'video', video, 'ipr', intraperiodDur, 'vidsize', '704x576', ...
        'alpha_q', alpha_q, 'alpha_f', alpha_f, 'fmax', 30, 'fr', [30 15], 'L', numLayers );    
    
    vs.maxValidTargetBitrate = 1300;
    vs.eta = zeros( size(vs.fr) );
    vs.theta = zeros( size(vs.fr) );
    vs = frameSizeModel( vs );
    
    % Derive qmin from the IPPP coding, we'll need this for comparison
    % qmin (achieving Q=1) is the QS at Rmax at 30 Hz with IPPP
    
    vs.L = 1;
    [actBitrate, q] = getQSandR( vs );
    [beta_q, ~, R0, q0] = getRSTARparam( actBitrate, q ) ;
    qmin = q0 * ((Rmax/R0))^(-1/beta_q);
    
    % Derive R-STAR parameters for numLayers temporal layers (hPP)
    vs.L = numLayers;
    [actBitrate, q] = getQSandR( vs );
    [beta_q, beta_f, R0, q0] = getRSTARparam( actBitrate, q );

    % at which rate do we have QS=qmin at 30 Hz? That's Rmax. (not really needed)
    vs.Rmax = R0 * (q0/qmin)^beta_q;
    
    vs.qmin = qmin;
    vs.beta_q = beta_q; vs.beta_f = beta_f;
    vs.R0 = R0; vs.q0 = q0;    
    
    fprintf('Sequence selected: %s\nR-STAR param: R0 = %f, beta_q = %f, q0 = %f, beta_f = %f, f0 = %f\n',...
        vs.video, vs.R0, vs.beta_q, vs.q0, vs.beta_f, 30);
    
    fprintf('Q-STAR param: qmin = %f, fmax = %f, Rmax = %f\n', vs.qmin, vs.fmax, vs.Rmax );
    fprintf('For R = %d kbps and f = 30 Hz, q = %f\n', Rmax, vs.q0 .* ( (Rmax./vs.R0) ).^(-1/vs.beta_q))

return

% TODO: QS calculation from QP is wrong
function [actBitrate, qs] = getQSandR( vs )
% getQSandR         parse files to get video bitrates and mean QS values
%  This function returns the encoded video bitrates and the mean QS for
%  each target bitrate and encoding frame rate for the sequence selected.
    
    fr = vs.fr; fr = fr(:);

    actBitrate = cell( size(fr) );
    qs = cell( size(fr) );
    
    for i = 1:length(fr)
        
        f = fr(i);
        
        datapath = strcat('data/',vs.video,'-',vs.vidsize,'-', ...
            strcat(num2str(f),'.0'),'-',num2str(vs.ipr*f),'-',num2str(vs.L),'/');
        
        [frameSizes, QPs, targetBitRates] = parse_log_files( datapath );    

        actBitrate{i} = sum(frameSizes)*0.008/10; % vid duration = 10 sec
        qs{i}         = mean( 2.^((QPs-4)./6) );
    end

    actBitrate = cell2mat(actBitrate);
    qs         = cell2mat(qs);
    
    % choose valid bitrates if need be
    valid = targetBitRates <= vs.maxValidTargetBitrate;
    
    actBitrate = actBitrate(:,valid);
    qs         = qs(:,valid);

return

function [beta_q, beta_f, R0, q0] = getRSTARparam( measuredBitrate, q_mean ) 
% Get RSTAR parameters with this function.
% TODO: code here assumes 2 temporal and 1 spatial resolutions, generalize
% TODO: return RSTARparam struct
    R0 = max(measuredBitrate(1,:));
    q0 = min(q_mean(1,:));
    
    ft = fittype('rstar( x, a, b )');
    model = fit((q0./q_mean(2,:))', (measuredBitrate(2,:)/R0)', ft, 'Startpoint', [0 0]);
    
    beta_q = model.a;
    beta_f = -log2(model.b);

return

function vs = frameSizeModel( vs )
    
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

        vs.theta( f == vs.fr ) = model.a;
        vs.eta( f == vs.fr ) = model.b;  
        
        fprintf('%s@%d Hz => theta=%f, eta=%f\n',vs.video,  f, model.a,  model.b);
    
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