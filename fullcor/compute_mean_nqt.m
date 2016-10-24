function [mean_quality, pmf] = compute_mean_nqt( param, markov, k, m )

pmf = compute_pmf( param.ref, k, m, precompute( markov, k, m ), markov.ss, 1, 0 );

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

function pmf = compute_pmf( ref, k, m, X, Pi, i, n )

N = length(k);
% recursion stopping condition
if i > N
    pmf = zeros(N+1,1);
    pmf(n+1) = pmf(n+1)+sum(Pi);
    return
end

% Case of arrv
pmf = compute_pmf( ref, k, m, X, X(i).T * X(i).Arrv * Pi, i+1, n+1 );

% Case of loss
% with the loss of this frame, we reach the end so don't factor in the jump prob
if ref(i)+i == N 
    pmf = pmf + compute_pmf( ref, k, m, X, X(i).Loss * Pi, N+1, n );
else
    ch = i+(1:ref(i));
    pmf = pmf + compute_pmf( ref, k, m, X, ...
        (X(i).T)^(1+sum(k(ch)+m(ch))) * X(i).Loss * Pi, i+1+ref(i), n );
end

return