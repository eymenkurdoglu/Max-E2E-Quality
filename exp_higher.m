function v = exp_higher(k)
G = length(k)/4;
k = (reshape(k,4,G))';
v = ones(G,1) + k(:,2) + k(:,3) + (k(:,3).*k(:,4));
return