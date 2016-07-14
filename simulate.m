function [fec_allocs, max_qual, max_fps] = simulate( alphas, betas, fec, l, qp, T, templayers )

N = size( l,1 ); % intra-period length
fec_allocs = zeros(N,length(fec),length(alphas));
max_qual = zeros(length(alphas),length(fec));
max_fps = zeros(length(alphas),length(fec));

for j = 1:length(alphas)

    alpha = alphas(j);
    beta  = betas(j);
    avg_num_frames  = zeros(1,length(fec));

    fprintf('\nSimulating for good-to-bad prob = %f\n', alpha);

    % search!
    for i = 1:length(fec)
        tic
        [best_alloc, best_obj] = greedy_fec_search3( fec(i), l(:,i), alpha, beta, templayers );
        toc
        avg_num_frames(i) = best_obj;
        fec_allocs(:,i,j) = best_alloc;
    end

    max_qual(j,:) = mnqq(mean(qp),4.51,20) .* mnqt(avg_num_frames/T,3.09,30);
    max_fps(j,:) = avg_num_frames/T;
end

return