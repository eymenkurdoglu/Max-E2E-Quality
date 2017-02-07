function [bestAlloc, bestNqt, bestPMF] = gfsmar( M, k, Markov, vs, bestAlloc )

k = k(:);

maxMeanNumDecFr = calcMeanNumDecFr( calcFrameArrvProbs( Markov, vs.referenceOf, k, bestAlloc ), vs.L, 1 );

while M > 0
    
    checkedFramesFrom = cell( 1, vs.L ); % don't check equal-sized frames from same layer
    
    for j = 1 : length(k)
        
        m = bestAlloc;
        
        thisLayer = vs.layerOf(j);
        if ~isempty( checkedFramesFrom{ thisLayer } ) && ... % < prevent 'exceeds mtx dimensions'
                any( checkedFramesFrom{ thisLayer }(1,:) == k(j) & checkedFramesFrom{ thisLayer }(2,:) == m(j) )
            continue;
        else
            checkedFramesFrom{ thisLayer } = [ checkedFramesFrom{ thisLayer }, [k(j);m(j)] ];
            m(j) = m(j) + 1;
        end

        value = calcMeanNumDecFr( calcFrameArrvProbs( Markov, vs.referenceOf, k, m ), vs.L, 1 );
        if value >= maxMeanNumDecFr
            maxMeanNumDecFr = value;
            pickedFrame = j;
        end
    end
    
    bestAlloc(pickedFrame) = bestAlloc(pickedFrame) + 1;
    
    M = M-1;
end

bestPMF = compPmf( vs.numDescendants, k, bestAlloc, precompute( Markov, k, bestAlloc ), Markov.ss, 1, 0 );

bestNqt = mnqt( (0:length(k))/vs.ipr, vs.alpha_f, vs.fmax ) * bestPMF;
return

function meanNumDecFr = calcMeanNumDecFr( p, L, topLevel )
%calcFrameRate   This function calculates the expected number of decoded
%frames, given the frame arrival probabilities and the number of layers in
%the hierarchical-P structure. Always call with topLevel==1 at highest
%level of recursion.

meanNumDecFr = 0;

if L == 0
    meanNumDecFr = 1;
else
    gopSize = 2^(L-1); % gop length
    
    numGops = length(p)/gopSize; % num of gops
    
    s = ( reshape(p, gopSize, numGops) )';
    
    for i = numGops:-1:1
        if i > 1
            meanNumDecFr = (meanNumDecFr + calcMeanNumDecFr( s(i,:), L-1, 0 )) * s(i,1);
        else
            meanNumDecFr = (meanNumDecFr + calcMeanNumDecFr( s(i,:), L-1, 0 ));
        end
    end
end

if topLevel
    meanNumDecFr = meanNumDecFr * p(1);
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

function p = calcFrameArrvProbs( markov, referenceOf, k, m )
%calcFrameArrvProbs   This function calculates the probability of arrival
%for each frame in the intra-period for the Markovian channels. The inputs
%are the frame sizes (k), number of FEC packets (m), the reference frame
%list for each frame, and finally the Markovian matrices. Output is a Nx1
%column vector of arrival probabilities.

p = zeros(2,length(k));

for i = 1:length(k)
    
    upperLimit = k(i)+m(i)-1;
    
    % what is Arrv?
    if upperLimit == 0 % <=> k(i)=1, m(i)=0 for this frame
        Arrv = [0,0;   ...
                0,1];
    else
        Arrv = [ sum( markov.L0(1 : m(i)-1, upperLimit) ), sum( markov.R0(k(i) : upperLimit, upperLimit) ); ...
                 sum( markov.L1(1 : m(i), upperLimit) ), sum( markov.R1(k(i)-1 : upperLimit, upperLimit) ) ];
    end
      
    refFrame = referenceOf(i);
    
    if refFrame == 0 % I-frame
        p(:,i) = Arrv * markov.ss;
    else
        inb = refFrame+1 : i-1; %inb(1) = []; inb(end) = []; !!!!!!!!!!!!!!!!!!!
        p(:,i) = Arrv * markov.T^(1 + sum( k(inb) + m(inb) )) * p(:,refFrame)/sum(p(:,refFrame));
    end
    
end
p = (sum(p))';
return