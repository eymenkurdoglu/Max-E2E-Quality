function p = calcFrameArrvProbs( markov, referenceMap, k, m )

N = length(k);
p = zeros(2,N);

for i = 1:N
    
    refFrame = referenceMap(i);
    upperLimit = k(i)+m(i)-1;
    
    if upperLimit == 0
        Arrv = [0,0;0,1];
    else
        Arrv = [ sum( markov.L0(1 : m(i)-1, upperLimit) ), sum( markov.R0(k(i) : upperLimit, upperLimit) ); 
        sum( markov.L1(1 : m(i), upperLimit) ), sum( markov.R1(k(i)-1 : upperLimit, upperLimit) ) ];
    end
      
    if refFrame < 0
        p(:,1) = Arrv * markov.ss;
    else
        inb = refFrame : i; 
        inb(1) = []; inb(end) = [];
        jumpOver = 1+sum( k(inb) + m(inb) );
        p(:,i) = Arrv * (markov.T)^jumpOver * p(:,refFrame)/sum(p(:,refFrame));
    end
    
end
p = (sum(p))';
return