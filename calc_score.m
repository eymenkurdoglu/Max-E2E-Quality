function score = calc_score( phi, mode )
% calc_score_hpp calculates the expected number of decodable frames.
% Success probabilities are directly taken from the "success_" vector.

% We will consider length(k) events (lost I-frame event is not considered
% as it has no effect on the expectation).

N = length(phi);

if mode == 0 % IPPP structure, this part calculates a + a*b + a*b*c + a*b*c*d + ...   
    score = 0;
    for i = N:-1:1
        score = (score+1)*phi(i);
    end
elseif mode == 1 % hPP structure
    G = N/4;
    s = ( reshape(phi, 4, G) )';
    b = s(:,1);
    e = ones(G,1) + s(:,2) + s(:,3) + (s(:,3).*s(:,4));
    
    score = 0;
    for i = G:-1:1
        score = (score+e(i))*b(i);
    end
end

return