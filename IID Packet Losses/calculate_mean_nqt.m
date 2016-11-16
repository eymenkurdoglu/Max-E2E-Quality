function [mean_quality, pmf] = calculate_mean_nqt( tree, phi, alpha_t, tmax, ipr )

pmf = calculate_pmf_ind( tree, phi );

mean_quality = mnqt( (0:length(pmf)-1)/ipr, alpha_t, tmax ) * pmf';

return

function pmf = calculate_pmf_ind( tree, phi )

node_index = tree.get(1);

if tree.isleaf(1)
    pmf = [ 1-phi( node_index ), phi( node_index ) ];
else
    a = [];
    for child = tree.getchildren(1)
        if isempty(a)
            a = calculate_pmf_ind( tree.subtree(child), phi );
        else
            a = conv( a, calculate_pmf_ind( tree.subtree(child), phi ) );
        end
    end
    pmf = [ 1-phi( node_index ), phi( node_index )*a ];
end

return