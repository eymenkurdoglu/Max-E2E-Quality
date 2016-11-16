function [mean_quality, pmf] = compute_mean_nqt( param, markov, k, m )

pmf = compute_pmf( param.numDescendants, k, m, precompute( markov, k, m ), markov.ss, 1, 0 );

mean_quality = mnqt( (0:length(k))/param.ipr, param.alpha_t, param.tmax ) * pmf;
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

function pmf = compute_pmf( numDescendants, k, m, X, Pi, i, n )

N = length(k);

if i > N % stop recursion
    pmf = zeros(N+1,1);
    pmf(n+1) = pmf(n+1)+sum(Pi);
    return
end

% Case of arrv
pmf = compute_pmf( numDescendants, k, m, X, X(i).T * X(i).Arrv * Pi, i+1, n+1 );

% Case of loss
if numDescendants(i) == N-i % rest of the frames undecodable
    pmf = pmf + compute_pmf( numDescendants, k, m, X, X(i).Loss * Pi, N+1, n );
else
    descendants = i+(1:numDescendants(i)); % could be empty
    jumpOver = 1+sum( k(descendants) + m(descendants) );
    pmf = pmf + compute_pmf( numDescendants, k, m, X, ...
        (X(i).T)^jumpOver * X(i).Loss * Pi, i+1+numDescendants(i), n );
end

return

function pmf = calculate_pmf( tree, k, m, markov, prob, i, N, children )

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
    
    pmf = calculate_pmf( tree, k, m, markov, A*B*prob, i+1, N, [] ) + ...
        calculate_pmf( tree.chop( node ), k, m, markov, C*B*prob, i+1, N, children );
else
    pmf = calculate_pmf( tree, k, m, markov, prob, i+1, N, children );
end

return