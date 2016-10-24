close all
clear all
clc
file = 'optimized.mat';
target = '~/Google Drive/NYU/Research/papers/fec/fig/';

if ~exist( file, 'file' )
    VIDEO = 'CREW';
    L = 3;
    P_GBs = [0.01, 0.001];
    P_BGs = 9*P_GBs;
    SRs  = 100 : 50 : 1300;
    
    optimalQ = zeros( length(P_GBs), length(SRs) );
    optimaleFR = zeros( length(P_GBs), length(SRs) );
    optimaldFR = zeros( length(P_GBs), length(SRs) );
    optimalRV = zeros( length(P_GBs), length(SRs) );
    optimalm = cell( length(P_GBs), length(SRs) );
    optimalPMF = cell( length(P_GBs), length(SRs) );
    
    for i = 1:length(P_GBs)
        p_gb = P_GBs(i);
        p_bg = P_GBs(i);
        fprintf('$$$ Simulating for P_GB = %f, P_BG = %f\n', p_gb,p_bg);

        for j = 1:length(SRs)
            SR = SRs(j);
            
            fprintf('$$$ - Sending rate = %d kbps\n', SR);
            
            [Qmax, FRopt, RVopt, mopt, PMFopt] = heuristic_search( p_gb, p_bg, SR, L, VIDEO );
            
            optimalQ(i,j) = Qmax;
            optimaleFR(i,j) = FRopt;
            optimalRV(i,j) = RVopt;
            optimalm{i,j} = mopt;
            optimalPMF{i,j} = PMFopt;
        end
    end

    save(file, 'optimalQ', 'optimaleFR', 'optimalRV', 'optimalm', 'optimalPMF', 'VIDEO', 'P_GBs', 'P_BGs', 'SRs', 'L')
else
    load( file )
end