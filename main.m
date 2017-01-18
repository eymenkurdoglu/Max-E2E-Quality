function main( sequences, pgb, pbg, bw )

intraperiodDur = 16/15;
numLayers = 3;

numChains = length(pgb);
numCapacs = length(bw);

if length(pbg) ~= numChains 
    display 'pbg and pgb lengths are not equal'
    return
elseif ~all(diff(bw) > 0)
    display 'Input bandwidths are not in right order'
    return
end

for video = sequences
    
    matFile = [video{1},'-',num2str(numLayers),'.mat'];

    if exist( matFile, 'file' )
        display 'Prior results found, move them somewhere else first'
        plotFigures( matFile )
        return
    end
        
    NQQ = zeros( numChains, numCapacs );
    NQT = zeros( numChains, numCapacs );
    F = zeros( numChains, numCapacs );
    R = zeros( numChains, numCapacs );
    M = cell ( numChains, numCapacs );
    D = cell ( numChains, numCapacs );
    
    vs = initstate(video,numLayers,bw(end),intraperiodDur);
    
    % plot the Q(R) for lossless case
    figure; hold on;
    plot( bw, max( mnqq( vs.q0 .* ( (bw./vs.R0).*(15/vs.fmax)^-vs.beta_f ).^(-1/vs.beta_q), vs.alpha_q, vs.qmin )...
        .* mnqt( 15, vs.alpha_f, vs.fmax ), mnqq( vs.q0 .* ( (bw./vs.R0) ).^(-1/vs.beta_q), vs.alpha_q, vs.qmin ) ) )
    
    for i = 1:numChains 
        
        piv = 1;
        fr = [30 15];
        
        for j = 1:numCapacs

            [nqt,nqq,f,r,m,d] = solve( pgb(i), pbg(i), bw(numCapacs-j+1), vs, fr, piv );
 
            NQQ(i,numCapacs-j+1) = nqq;
            NQT(i,numCapacs-j+1) = nqt;
            F(i,numCapacs-j+1) = f;
            R(i,numCapacs-j+1) = r;
            M{i,numCapacs-j+1} = m;
            D{i,numCapacs-j+1} = d;
            
            piv = 1-((1-r/bw(numCapacs-j+1))*0.6);
            fprintf('### FEC perc = %f, pivot set to %f\n',100*(1-r/bw(numCapacs-j+1)),piv)
            
            if f < 30
                fr = 15;
            end
        end
    end
    save(matFile,'NQQ','NQT','F','R','M','D','pgb','pbg','bw','numLayers','vs')
end

return

function vs = initstate(video,numLayers,Rmax,intraperiodDur)
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
        'alpha_q', alpha_q, 'alpha_f', alpha_f, 'fmax', 30, 'fr', [30 15] );    
    
    vs.maxValidTargetBitrate = 1300;
    
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

function plotFigures( matFile )

load(matFile)

numCurves = length(pgb);
l = cell(1,numCurves); % legend strings
for i = 1:numCurves
    l{i}=['PLR=',num2str(pgb(i))];
end

figure
plot(bw,(NQQ.*NQT)');
xlabel('Sending Bitrate (Mbps)'); legend(l);

figure
plot(bw,100*(1-(R./repmat(bw,numCurves,1))')); legend(l);

figure
plot(bw,F'); legend(l);

return