function AverageLengths = mean_length_per_layer( Lengths, N, L )

AverageLengths = cell(1,L+1);

for l = L+1 : -1 : 3
    AverageLengths{l} = Lengths(2:2:end);
    Lengths(2:2:end) = [];
end

AverageLengths{1} = Lengths(1:(N*2^(1-L)):end);

Lengths(1:(N*2^(1-L)):end) = [];

AverageLengths{2} = Lengths;

AverageLengths = (cellfun(@mean, AverageLengths))';

return