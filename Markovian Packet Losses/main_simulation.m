close all
clc
VIDEOS = {'CREW','CITY','HARBOUR','FG'};

for v = 1:length(VIDEOS)
    
video = VIDEOS{v};
display(video);
    
L = 3;
file = ['optimized',video,'.mat'];

if ~exist( file, 'file' )

    P_GBs = [0.01, 0.001];
    P_BGs = 9*P_GBs;
    SRs  = 100 : 50 : 1300;
    
    optimalQ = zeros( length(P_GBs), length(SRs) );
    optimaleFR = zeros( length(P_GBs), length(SRs) );
    optimalRV = zeros( length(P_GBs), length(SRs) );
    optimalm = cell( length(P_GBs), length(SRs) );
    optimalPMF = cell( length(P_GBs), length(SRs) );
    
    for i = 1:length(P_GBs)
        p_gb = P_GBs(i);
        p_bg = P_BGs(i);
        
        fprintf('# Simulating for P_GB = %f, P_BG = %f\n', p_gb,p_bg);

        tryboth = 1;
        for j = 1:length(SRs)
            sendingRate = SRs(j);
            
            fprintf('######### Sending rate = %d kbps\n', sendingRate);
            
            [Qmax, FRopt, RVopt, mopt, PMFopt] = heuristic_search( p_gb, p_bg, sendingRate, L, video, tryboth );
            if FRopt == 30
                tryboth = 0;
            end
            optimalQ(i,j) = Qmax;
            optimaleFR(i,j) = FRopt;
            optimalRV(i,j) = RVopt;
            optimalm{i,j} = mopt;
            optimalPMF{i,j} = PMFopt;
        end
    end

    save(file, 'optimalQ', 'optimaleFR', 'optimalRV', 'optimalm', 'optimalPMF', 'video', 'P_GBs', 'P_BGs', 'SRs', 'L')
else
    load( file )
end

end