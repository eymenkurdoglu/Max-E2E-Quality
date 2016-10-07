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