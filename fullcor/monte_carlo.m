function nqt = monte_carlo( markov, k, m, tree, param, n )
N = length(k);
len = sum(k+m);
nqt = 0;
T = tree;
alpha = markov.alpha;
beta = markov.beta;
for sample = 1:n
    bitstr = ones(len,1);
    if rand < alpha/(alpha+beta)
        bitstr(1) = 0;
    end

    if len > 1
        for i = 2:len
            r = rand;
            if ( bitstr(i-1) == 1 && r < alpha ) || ( bitstr(i-1) == 0 && r < 1-beta )
                bitstr(i) = 0;
            end
        end
    end

    t = T;
    end_ = 0;
    for i = 1:N
        start_ = end_+1;
        end_ = start_ + k(i) + m(i) - 1;
        if any( t == i )
            if sum(bitstr(start_:end_)>0) < k(i)
                t = t.chop( find(t==i) );
            end
        end
    end

    nqt = nqt + mnqt(t.nnodes/param.ipr,param.alpha_t,param.tmax);
end
nqt = nqt/n;
return