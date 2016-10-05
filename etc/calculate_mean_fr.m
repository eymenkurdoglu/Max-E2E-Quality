function score = calculate_mean_fr( phi, L, top )

score = 0;

if L == 0
    score = 1;
else
    g = 2^(L-1); % gop length
    
    G = length(phi)/g; % num of gops
    
    s = ( reshape(phi, g, G) )';
    
    for i = G:-1:1
        if i > 1
            score = (score + calculate_mean_fr( s(i,:), L-1, 0 )) * s(i,1);
        else
            score = (score + calculate_mean_fr( s(i,:), L-1, 0 ));
        end
    end 
end

if top
    score = score * phi(1);
end

return