function [bestAlloc, maxValue, bestPMF] = greedy_fec_search( M, k, Markov, Param, bestAlloc )

k = k(:);

% [maxValue, bestPMF] = compute_mean_nqt( Param, Markov, k, bestAlloc );
% best_nqt = monte_carlo( markov, k, best_alloc, tree, param, 2000 );
maxValue = calcFrameRate( calcFrameArrvProbs( Markov, Param.referenceMap, k, bestAlloc ), Param.numTempLayers, 1 );
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

%         nqt = monte_carlo( markov, k, m, tree, param, 2000 );
%         [value, pmf] = compute_mean_nqt( Param, Markov, k, m ); % precomputes matrices inside
        value = calcFrameRate( calcFrameArrvProbs( Markov, Param.referenceMap, k, m ), Param.numTempLayers, 1 );
        if value >= maxValue
%             bestPMF = pmf;
            maxValue = value;
            pickedFrame = j;
        end
    end
    
    bestAlloc(pickedFrame) = bestAlloc(pickedFrame) + 1;
    
    M = M-1;
end
[maxValue, bestPMF] = compute_mean_nqt( Param, Markov, k, bestAlloc );
return