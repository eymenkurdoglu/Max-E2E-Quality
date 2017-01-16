function Q = mnqt( t, alpha_t, tmax )

num = 1-exp(-1*alpha_t*(t./tmax).^0.63);
denom = 1-exp(-1*alpha_t);

Q = num/denom;

return