function score = calc_score_x( phi, layer, top )

score = 0;

if layer == 0
    
    score = 1;
    
else
    g = 2^(layer-1); % gop length
    
    G = length(phi)/g; % num of gops
    
    s = ( reshape(phi, g, G) )';
    
    for i = G:-1:1
        if i > 1
            score = (score + calc_score_x( s(i,:), layer-1, 0 )) * s(i,1);
        else
            score = (score + calc_score_x( s(i,:), layer-1, 0 ));
        end
    end 
end

if top
    score = score * phi(1);
end

return