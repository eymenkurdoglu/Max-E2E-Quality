close all
clear all
file = 'optimized.mat';
target = '~/Google Drive/NYU/Research/papers/fec/fig/';

if ~exist( file, 'file' )
    VIDEO = 'CREW';
    L = 3;
    PLRs = 0.05 : 0.05 : 0.25;
    SRs  = 100 : 50 : 1300;
    
    optimalQ = zeros( length(PLRs), length(SRs) );
    optimaleFR = zeros( length(PLRs), length(SRs) );
    optimaldFR = zeros( length(PLRs), length(SRs) );
    optimalRV = zeros( length(PLRs), length(SRs) );
    optimalm = cell( length(PLRs), length(SRs) );
    
    for i = 1:length(PLRs)
        PLR = PLRs(i);
        
        fprintf('# Simulating for PLR = %f\n', PLR);

        for j = 1:length(SRs)
            SR = SRs(j);
            
            fprintf('######### Sending rate = %d kbps\n', SR);
            
            [Qmax, FRopt, RVopt, mopt] = heuristic_search( PLR, SR, L, VIDEO );
            
            optimalQ(i,j) = Qmax;
            optimaleFR(i,j) = FRopt;
            optimalRV(i,j) = RVopt;
            optimalm{i,j} = mopt;
        end
    end

    save(file, 'optimalQ', 'optimaleFR', 'optimalRV', 'optimalm', 'VIDEO', 'PLRs', 'SRs', 'L')
else
    load( file )
end