function y = fitmodel( x, a, b )

y = 1-exp(-1*a*x.^b);

return