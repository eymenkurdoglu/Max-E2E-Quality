function main( video, pgb, pbg, bw, numLayers )

assert( length(pbg) == length(pgb) )
assert( all(diff(bw) > 0) )

PACKET_SIZE = 200;
INTRAPERIOD = 16/15;

numChains = length(pgb);
numCapacs = length(bw);
  
matFile = ['results/',video,'-',num2str(numLayers),'-15Hz.mat'];

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

vs = initState( video, numLayers, bw(end), INTRAPERIOD );

dbstop if error

parfor i = 1:numChains 

    piv = 1; fr = [30 15]; BW = bw;
    
    % auxiliary variables for parfor 
    NQQ_ = zeros( 1, numCapacs );
    NQT_ = zeros( 1, numCapacs );
    F_ = zeros( 1, numCapacs );
    R_ = zeros( 1, numCapacs );
    M_ = cell ( 1, numCapacs );
    D_ = cell ( 1, numCapacs );
    
    for j = numCapacs : -1 : 1 % start from the greatest bw

        [nqt,nqq,f,r,m,d] = solve( pgb(i), pbg(i), BW(j), vs, fr, piv );

        NQQ_(j) = nqq; NQT_(j) = nqt; F_(j) = f; R_(j) = r; M_{j} = m; D_{j} = d;

        if f == 15; fr = 15; end
        piv = 1 - (1-r/bw(j))*0.6;
        
%         fprintf('### FEC perc = %f, pivot set to %f\n',100*(1-r/bw(j)),piv)
    end
    
    NQQ(i,:) = NQQ_;
    NQT(i,:) = NQT_;
    F(i,:) = F_;
    R(i,:) = R_;
    M(i,:) = M_;
    D(i,:) = D_;
    
%     fprintf( '> Done with pgb = %f...\n', pgb(i) );
end

save(matFile,'NQQ','NQT','F','R','M','D','pgb','pbg','bw','numLayers','vs','PACKET_SIZE')

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