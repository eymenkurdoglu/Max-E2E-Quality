function [best_alloc, best_nqt, best_pmf] = gfsiid( M, k, PLR, vs, varargin )

k = k(:);
N = length(k);

if ~isempty(varargin)
    in = varargin{1}; in = in(:);
else
    in = zeros(length(k),1);
end

if PLR == 0
    best_alloc = in;
    best_nqt = N;
    best_pmf = zeros(1,N+1); best_pmf(end) = 1;
    return
end

p = zeros(N,2);
for i = 1:N
    % First column : frame arrival probabilities
    p(i,1) = sum(binopdf( 0:in(i)  , in(i)+k(i)  , PLR));
    % Second column: frame arrival probabilities with an additional FEC block
    p(i,2) = sum(binopdf( 0:in(i)+1, in(i)+k(i)+1, PLR));
end

tree = vs.t;

best_alloc = in;
[best_nqt, best_pmf] = calculate_mean_nqt( tree, p(:,1), vs.alpha_f, vs.fmax, vs.ipr );

while M > 0
    % check each neighbor
    for j = randperm(N)
        
        p_    = p(:,1);
        p_(j) = p(j,2);
        
        [nqt, pmf] = calculate_mean_nqt( tree, p_, vs.alpha_f, vs.fmax, vs.ipr );
        
        if nqt >= best_nqt
            best_pmf = pmf;
            best_nqt = nqt;
            picked = j;
        end
    end
    
    best_alloc(picked) = best_alloc(picked)+1;
    
    % need to only change the success prob corr. to the picked frame
    p(picked,1) = p(picked,2);
    
    % to update the value in the second column
    p(picked,2) = sum( binopdf(0:(best_alloc(picked)+1), k(picked)+best_alloc(picked)+1, PLR) );
    
    M = M-1;
end

return

function [mean_quality, pmf] = calculate_mean_nqt( tree, p, alpha_t, fmax, ipr )

pmf = calculate_pmf_ind( tree, p );

mean_quality = mnqt( (0:length(pmf)-1)/ipr, alpha_t, fmax ) * pmf';

return

function pmf = calculate_pmf_ind( t, p )

node_index = t.get(1);

if t.isleaf(1)
    pmf = [ 1-p( node_index ), p( node_index ) ];
else
    a = [];
    for child = t.getchildren(1)
        if isempty(a)
            a = calculate_pmf_ind( t.subtree(child), p );
        else
            a = conv( a, calculate_pmf_ind( t.subtree(child), p ) );
        end
    end
    pmf = [ 1-p( node_index ), p( node_index )*a ];
end

return