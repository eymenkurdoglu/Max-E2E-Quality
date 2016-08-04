dbstop if error
clc
eps = 0.25;
k = [15 14];
maks= 10;

M = 0:maks;
y = zeros(1,length(M));
for m = M;
    y(m+1) = sum( binopdf(0:m, k(1)+m, eps) ) + sum( binopdf(0:(maks-m), k(2)+maks-m, eps) );
end