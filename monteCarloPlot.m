function [MeanFrameIntervals, StdFrameIntervals] = monteCarloPlot(video,L,numRun)

dbstop if error

load( [video,'-',num2str(L),'.mat'] )

numChains = length(pgb);
numCapacs = length(bw);

PGB = pgb;
PBG = pbg;

T(1) =  create_hierP_tree( L, 16 );
T(2) =  create_hierP_tree( L, 32 );

MeanFrameIntervals = zeros(numChains,numCapacs);
StdFrameIntervals = zeros(numChains,numCapacs);

K = cell(numChains,numCapacs);
for i = 1:numChains 
    for j = 1:numCapacs
        vs.f = F(i,j);
        k = ceil( estimFrameSz( vs, R(i,j) ) / PACKET_SIZE );
        K{i,j} = k;
    end
end

parfor i = 1:numChains 
    
    p_gb = PGB(i); p_bg = PBG(i); iid = 0;
    if p_gb + p_bg == 1; iid = 1; end
    
    for j = 1:numCapacs
        
        k = K{i,j};
        m = M{i,j};    
        
        assert( length(m) == length(k) ); N = length(m);
        
        meanFrameInterval = 0; stdFrameInterval = 0;
        
        for t = 1:numRun % mc
            
            if iid; packetLossPattern = rand(1,sum(k+m)) > p_gb; end
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
            intervals = diff( sort( [cell2mat( tree.Node ); N+1] ) );
            if ~isempty( intervals ); meanFrameInterval = meanFrameInterval + mean(intervals); 
                stdFrameInterval = stdFrameInterval + std(intervals); end
        end
        MeanFrameIntervals(i,j) = meanFrameInterval / numRun;
        StdFrameIntervals(i,j) = stdFrameInterval / numRun;
    end
end

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