function Q = mnqq( q, alpha_q, qmin )

num = 1-exp(-1*alpha_q*(qmin./q));
denom = 1-exp(-1*alpha_q);

Q = num/denom;

return