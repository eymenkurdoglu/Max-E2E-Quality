function framerate = calcFrameRate( p, L, topLevelFlag )

framerate = 0;

if L == 0
    framerate = 1;
else
    g = 2^(L-1); % gop length
    
    G = length(p)/g; % num of gops
    
    s = ( reshape(p, g, G) )';
    
    for i = G:-1:1
        if i > 1
            framerate = (framerate + calcFrameRate( s(i,:), L-1, 0 )) * s(i,1);
        else
            framerate = (framerate + calcFrameRate( s(i,:), L-1, 0 ));
        end
    end 
end

if topLevelFlag
    framerate = framerate * p(1);
end

return