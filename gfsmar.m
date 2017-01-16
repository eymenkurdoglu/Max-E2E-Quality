function [bestAlloc, bestNqt, bestPMF] = gfsmar( M, k, Markov, Param, bestAlloc )

k = k(:);

maxFrameRate = calcFrameRate( calcFrameArrvProbs( Markov, Param.referenceMap, k, bestAlloc ), ...
    Param.numTempLayers, 1 );

while M > 0
    
    checkedFramesFrom = cell( 1, Param.numTempLayers ); % don't check equal-sized frames from same layer
    
    for j = 1 : length(k)
        
        m = bestAlloc;
        
        thisLayer = Param.layerMap(j);
        if ~isempty( checkedFramesFrom{ thisLayer } ) && ... % < prevent 'exceeds mtx dimensions'
                any( checkedFramesFrom{ thisLayer }(1,:) == k(j) & checkedFramesFrom{ thisLayer }(2,:) == m(j) )
            continue;
        else
            checkedFramesFrom{ thisLayer } = [ checkedFramesFrom{ thisLayer }, [k(j);m(j)] ];
            m(j) = m(j) + 1;
        end

        value = calcFrameRate( calcFrameArrvProbs( Markov, Param.referenceMap, k, m ), Param.numTempLayers, 1 );
        if value >= maxFrameRate
            maxFrameRate = value;
            pickedFrame = j;
        end
    end
    
    bestAlloc(pickedFrame) = bestAlloc(pickedFrame) + 1;
    
    M = M-1;
end

bestPMF = compPmf( Param.numDescendants, k, bestAlloc, precompute( Markov, k, bestAlloc ), Markov.ss, 1, 0 );

bestNqt = mnqt( (0:length(k))/Param.ipr, Param.alpha_t, Param.tmax ) * pmf;
return

function framerate = calcFrameRate( p, L, topLevel )
% always call with topLevel==1 at highest level of recursion

framerate = 0;

if L == 0
    framerate = 1;
else
    g = 2^(L-1); % gop length
    
    G = length(p)/g; % num of gops
    
    s = ( reshape(p, g, G) )';
    
    for i = G:-1:1
        if i > 1
            framerate = (framerate + calcFrameRate( s(i,:), L-1, 0 )) * s(i,1);
        else
            framerate = (framerate + calcFrameRate( s(i,:), L-1, 0 ));
        end
    end 
end

if topLevel
    framerate = framerate * p(1);
end

return

function X = precompute( markov, k, m )

N = length(k);
Arrv = cell(1,N);
Loss = cell(1,N);

for i = 1:N
    ul = k(i)+m(i)-1;
    if ul == 0
        Arrv{i} = [0,0;0,1]; Loss{i} = [1,0;0,0];
    else
    Arrv{i} = [ sum( markov.L0(1 : m(i)-1, ul) ), sum( markov.R0(k(i) : ul, ul) ); 
        sum( markov.L1(1 : m(i), ul) ), sum( markov.R1(k(i)-1 : ul, ul) ) ];
    Loss{i} = [ sum( markov.L0(max(m(i),1) : ul, ul) ), sum( markov.R0(1 : k(i)-1, ul) );
        sum( markov.L1(m(i)+1 : ul, ul) ), sum( markov.R1(1 : k(i)-2, ul) ) ];
    end
end

X = struct('Arrv',Arrv,'Loss',Loss,'T',markov.T);

return

function pmf = compPmf( numDescendants, k, m, X, Pi, i, n )

N = length(k);

if i > N % stop recursion
    pmf = zeros(N+1,1);
    pmf(n+1) = pmf(n+1)+sum(Pi);
    return
end

% Case of arrv
pmf = compPmf( numDescendants, k, m, X, X(i).T * X(i).Arrv * Pi, i+1, n+1 );

% Case of loss
if numDescendants(i) == N-i % rest of the frames undecodable
    pmf = pmf + compPmf( numDescendants, k, m, X, X(i).Loss * Pi, N+1, n );
else
    descendants = i+(1:numDescendants(i)); % could be empty
    jumpOver = 1+sum( k(descendants) + m(descendants) );
    pmf = pmf + compPmf( numDescendants, k, m, X, ...
        (X(i).T)^jumpOver * X(i).Loss * Pi, i+1+numDescendants(i), n );
end

return

function pmf = compPmfTree( tree, k, m, markov, prob, i, N, children )

% recursion stopping condition: reaching the end of the intra-period
if i > N
    pmf = zeros(N+1,1);
    pmf( tree.nnodes+1 ) = pmf( tree.nnodes+1 ) + sum(prob);
    return
end

if any( tree == i )
    ul = k(i)+m(i)-1;
    node = find( tree == i );
    
    A = [ sum( markov.L0(1 : m(i)-1, ul) ), sum( markov.R0(k(i) : ul, ul) ); 
        sum( markov.L1(1 : m(i), ul) ), sum( markov.R1(k(i)-1 : ul, ul) ) ];
    B = markov.T^(sum( k(children)+m(children) )+1);
    C = [ sum( markov.L0(max(m(i),1) : ul, ul) ), sum( markov.R0(1 : k(i)-1, ul) );
        sum( markov.L1(m(i)+1 : ul, ul) ), sum( markov.R1(1 : k(i)-2, ul) ) ];
    
    children = cell2mat( tree.subtree( node ).Node ); children(1) = [];
    
    pmf = compPmfTree( tree, k, m, markov, A*B*prob, i+1, N, [] ) + ...
        compPmfTree( tree.chop( node ), k, m, markov, C*B*prob, i+1, N, children );
else
    pmf = compPmfTree( tree, k, m, markov, prob, i+1, N, children );
end

return