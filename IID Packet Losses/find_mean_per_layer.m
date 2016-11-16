function mean_x = find_mean_per_layer( x, N, L )

mean_x = cell(1,L+1);

for l = L+1 : -1 : 3
    mean_x{l} = x(2:2:end);
    x(2:2:end) = [];
end

mean_x{1} = x(1:(N*2^(1-L)):end);

x(1:(N*2^(1-L)):end) = [];

mean_x{2} = x;

mean_x = (cellfun(@mean, mean_x))';

return