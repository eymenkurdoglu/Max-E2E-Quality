function [NQTmax, NQQmax, fopt, Ropt, mopt, PMFopt] = solve( p_gb, p_bg, bw, vs, fr, piv )

PACKET_SIZE = 200;
SEARCH_PERSISTENCE = 20;

% fprintf('##### Sending rate = %d kbps\n', bw)

% total #packets at this sending rate
totNumPack = floor( bw * vs.ipr / (0.008*PACKET_SIZE) );

%% Start the search
Qmax = 0;

for f = fr
    
    vs.f = f;
    vs.N = vs.f * vs.ipr;
    
    vs = intraPeriodStruct( vs );
    
    if p_bg + p_gb ~= 1
        markovFileName = ['markov-',num2str(1/p_bg),'-',num2str(p_gb/(p_gb+p_bg)),'.mat'];
        if exist( markovFileName, 'file' )
            load(markovFileName)
            assert( size(markov.L0,1) >= totNumPack )
        else
            markov = create_prob_struct( p_gb, p_bg, totNumPack );
            save(markovFileName, 'markov')
        end
    end
    
    prevk = zeros(vs.N,1); descending = 0;
    m = zeros(vs.N,1); 
    
    highest_Q_so_far = 0; 
    best_rate_so_far = 0; 
    best_pmf_so_far = 0;
    highest_NQT_so_far = 0; 
    highest_NQQ_so_far = 0;
    
    for R = ceil(bw * piv) : -1 : 50
        
        k = ceil( estimFrameSz(vs,R) / PACKET_SIZE );
        
        if all( k == prevk ) || totNumPack < sum(k)
            continue; 
        else
            prevk = k;
        end
        
%         fprintf('%d Hz, %d kbps: ', vs.f, R)
        
        if p_bg + p_gb ~= 1
            [m, NQT, pmf] = gfsmar( totNumPack-sum(k+m), k, markov, vs, m );
        else
            [m, NQT, pmf] = gfsiid( totNumPack-sum(k+m), k, p_gb, vs, m );
        end
        
        q = vs.q0 * ( (R/vs.R0)*(vs.f/vs.fmax)^-vs.beta_f )^(-1/vs.beta_q);
        NQQ = mnqq(q,vs.alpha_q,vs.qmin);     
        Q = NQQ * NQT;
        
%         fprintf('NQT = %f, NQQ = %f, Q = %f ',NQT,NQQ,Q);
        
        if Q > highest_Q_so_far
%             fprintf('++\n');
            descending = 0;
            best_pmf_so_far = pmf;
            highest_Q_so_far = Q;
            highest_NQT_so_far = NQT;
            highest_NQQ_so_far = NQQ;
            best_rate_so_far = R;
            best_alloc_so_far = m;
        else
%             fprintf('--\n');
            descending = descending + 1;
            if descending > SEARCH_PERSISTENCE
%                 fprintf('>> R = %d selected for f = %d\n', best_rate_so_far, f);
                break;
            end
        end
    end

    if highest_Q_so_far > Qmax
        Qmax = highest_Q_so_far;
        NQTmax = highest_NQT_so_far;
        NQQmax = highest_NQQ_so_far;
        PMFopt = best_pmf_so_far;
        fopt = f;
        Ropt = best_rate_so_far;
        mopt = best_alloc_so_far;
    end
%     fprintf('### %d Hz selected and R = %d\n', fopt, Ropt);
end

return

function vs = intraPeriodStruct( vs )

    L = vs.L;
    N = vs.N;
    
    numDescendants = zeros(1,N); % number of descendants
    gop = 2^(L-1);
    all_indices_so_far = 1 : gop : N; % base layer frames
    
    for i = all_indices_so_far
        numDescendants(i) = N-i;
    end
    
    while length(all_indices_so_far) < N
        gop = gop/2;
        numDescendants( all_indices_so_far + gop ) = gop - 1;
        all_indices_so_far = [all_indices_so_far, all_indices_so_far + gop];
    end

    vs.numDescendants = numDescendants;

    layerOf = zeros(N,1);
    for layer = L : -1 : 1   
        layerOf(mod(0:N-1,2^(L-layer))==0) = layer;
    end

    vs.layerOf = layerOf;
    vs.t = create_hierP_tree( L, N );
    vs.referenceOf = createReferenceMap(vs.t);
    
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

function referenceMap = createReferenceMap(tree)

referenceMap = zeros(tree.nnodes,1);

for i = 1:tree.nnodes
    d = tree.Parent(find(tree==i));
    if d == 0
        assert(i==1);
        referenceMap(i) = 0;
    else
        referenceMap(i) = tree.Node{d};
    end
end

return

function P = create_prob_struct( p_gb, p_bg, n )

% "0" and "1" duration tail distribution functions and probability mass functions
P_R = [1,p_bg*(1-p_gb).^(0:n-2)];
P_L = [1,p_gb*(1-p_bg).^(0:n-2)];
p_R = [1-p_bg,p_gb*P_R(2:end)];
p_L = [1-p_gb,p_bg*P_L(2:end)];

L  = create_mtx( n, P_R, p_R ); % Pr(m-1 losses occur in the next n-1 following a loss)
L0 = create_mtx( n, p_R, p_R ); % Pr(m-1 losses occur in the next n-1 between 2 losses)
L1 = L-L0; % Pr(m-1 losses occur in the next n-1 between a loss and a reception)

R  = create_mtx( n, P_L, p_L ); % Pr(m-1 receptions occur in the next n-1 following a reception)
R1 = create_mtx( n, p_L, p_L ); % Pr(m-1 receptions occur in the next n-1 between 2 receptions)
R0 = R-R1; % Pr(m-1 receptions occur in the next n-1 between a reception and a loss)

T = [1-p_bg,   p_gb;
       p_bg, 1-p_gb];
   
ss = [p_gb;p_bg]/(p_gb+p_bg);

P = struct('L0',L0,'L1',L1,'R0',R0,'R1',R1,'T',T,'ss',ss,'alpha',p_gb,'beta',p_bg);

return

function A = create_mtx( n, B, b )
% Make sure B and b have correct length (row vectors of size(A,1))
% B: row vector st B(i)=Pr(0^(i-1) | 1) or Pr(1^(i-1) | 0)
% b: row vector st b(i)=Pr(0^(i-1)1 | 1) or Pr(1^(i-1)0 | 0)
A = zeros(n);
A(1,:) = B;
for i = 2:n
    for j = i:n
        A(i,j) = fliplr(b(1:(j-i+1))) * A(i-1,i-1:j-1)';
    end
end

return