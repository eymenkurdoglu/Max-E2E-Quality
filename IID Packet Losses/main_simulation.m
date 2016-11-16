close all
clear all
clc
VIDEOS = {'CREW','CITY','HARBOUR','FG'};

for v = 1:length(VIDEOS)

video = VIDEOS{v};
display(video);
    
L = 3;
file = ['optimized',video,'.mat'];

if ~exist( file, 'file' )
    
    PLRs = 0.05 : 0.05 : 0.25;
    SRs  = 100 : 50 : 1300;
    
    optimalQ = zeros( length(PLRs), length(SRs) );
    optimaleFR = zeros( length(PLRs), length(SRs) );
    optimalRV = zeros( length(PLRs), length(SRs) );
    optimalm = cell( length(PLRs), length(SRs) );
    optimalPMF = cell( length(PLRs), length(SRs) );
    
    for i = 1:length(PLRs)
        PLR = PLRs(i);
        
        fprintf('# Simulating for PLR = %f\n', PLR);

        for j = 1:length(SRs)
            SR = SRs(j);
            
            fprintf('######### Sending rate = %d kbps\n', SR);
            
            [Qmax, FRopt, RVopt, mopt, PMFopt] = heuristic_search( PLR, SR, L, video );
            
            optimalQ(i,j) = Qmax;
            optimaleFR(i,j) = FRopt;
            optimalRV(i,j) = RVopt;
            optimalm{i,j} = mopt;
            optimalPMF{i,j} = PMFopt;
        end
    end

    save(file, 'optimalQ', 'optimaleFR', 'optimalRV', 'optimalm', 'optimalPMF', 'video', 'PLRs', 'SRs', 'L')
else
    load( file )
end

end