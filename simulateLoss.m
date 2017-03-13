function simulateLoss( video, L, numRun )

dbstop if error

isMarkovian = 1;

if isMarkovian; channel = '-markov'; else channel = ''; end
load( ['results/',video,'-',num2str(L),channel,'.mat'] )

numChains = length(pgb);
numCapacs = length(bw);

PGB = pgb;
PBG = pbg;
F_ = F;
M_ = M;

F = F_;
M = M_;

T(1) =  create_hierP_tree( L, 16 );
T(2) =  create_hierP_tree( L, 32 );

MeanFrameIntervals = zeros(numChains,numCapacs);
StdFrameIntervals = zeros(numChains,numCapacs);
MeanNumDecFrames = zeros(numChains,numCapacs);
numFreezes = zeros(numChains,numCapacs);

for i = 1:numChains 
    
    p_gb = PGB(i); p_bg = PBG(i);
    
    for j = 1:numCapacs
        
        k = K{i,j};
        m = M{i,j};    
        
        assert( length(m) == length(k) ); N = length(m);
        
        meanFrameInterval = 0; stdFrameInterval = 0; MeanNumDecFrame = 0; fullEmpty = 0;
        
        for t = 1:numRun % mc
            
            if ~isMarkovian; packetLossPattern = rand(1,sum(k+m)) > p_gb;
            else
                packetLossPattern = ones(1,sum(k+m));
                packetLossPattern(1) = rand > p_gb/(p_gb+p_bg);
                for p = 2 : sum(k+m)
                    r = rand;
                    if ( packetLossPattern(p-1) == 1 && r < p_gb ) || ( packetLossPattern(p-1) == 0 && r < 1-p_bg )
                        packetLossPattern(p) = 0;
                    end
                end
            end
            
            if F(i,j) == 15; tree = T(1); else tree = T(2); end
            
            end_ = 0;
            for u = 1:N
                start_ = end_+1;
                end_ = start_ + k(u) + m(u) - 1;
                if any( tree == u )
                    if sum(packetLossPattern(start_:end_)>0) < k(u)
                        tree = tree.chop( find(tree == u) );
                    end
                end
            end
            
            % gather stats
            intervals = diff( sort( [cell2mat( tree.Node ); N+1] ) ) / F(i,j);

            if ~isempty( intervals ); 
                pdf = intervals'/sum(intervals);
                timeAvgInterval = pdf * intervals;
                timeStdInterval = sqrt( pdf * ((intervals-timeAvgInterval).^2) );

                meanFrameInterval = meanFrameInterval + timeAvgInterval; 
                stdFrameInterval = stdFrameInterval + timeStdInterval; 
                MeanNumDecFrame = MeanNumDecFrame + tree.nnodes;
            else
                fullEmpty = fullEmpty + 1;
            end
        end
        MeanFrameIntervals(i,j) = meanFrameInterval / numRun;
        StdFrameIntervals(i,j) = stdFrameInterval / numRun;
        MeanNumDecFrames(i,j) = MeanNumDecFrame / numRun;
        numFreezes(i,j) = fullEmpty;
    end
end
save(['results/',video,'-',num2str(L),channel,'-MC_deneme.mat'],'MeanFrameIntervals','StdFrameIntervals','MeanNumDecFrames','numFreezes')
return

function t = create_hierP_tree( L, N )

gop = 2^(L-1);
G = N / gop;
t = tree(1);
for b = 2:G
    t = t.addnode(b-1,1+(b-1)*gop);
end

while nnodes(t) ~= N
    gop = gop / 2;
    for node = t.nodeorderiterator
        index = t.get( node );
        t = t.addnode( node, index + gop );
    end
end

return