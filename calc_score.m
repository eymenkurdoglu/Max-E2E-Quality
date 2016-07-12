function score = calc_score( success_, mode )
% calc_score_hpp calculates the expected number of decodable frames.
% Success probabilities are directly taken from the "success_" vector.

% We will consider length(k) events (lost I-frame event is not considered
% as it has no effect on the expectation).

N = length(success_);

if mode == 0 % IPPP structure, this part calculates a + a*b + a*b*c + a*b*c*d + ...
    score = 0;
    for i = N:-1:1
        score = (score+1)*success_(i);
    end
elseif mode == 1 % hPP structure
    gop = 4;
    assert( mod(N,gop)==0 );

    G = N/gop; % number of gops assuming gop length is 4
    base = fliplr(1:4:N);

    success_gop = (reshape(success_, 4, G))';
    v = ones(G,1) + success_gop(:,2) + success_gop(:,3) + (success_gop(:,3).*success_gop(:,4));
    cum_v = cumsum(v);

    score = cum_v(end);
    j = G-1;
    for i = base
        if j > 0
            score = cum_v(j) + (score-cum_v(j))*success_(i);
        else
            score = score*success_(i);
        end
        j = j-1;
    end
end

return