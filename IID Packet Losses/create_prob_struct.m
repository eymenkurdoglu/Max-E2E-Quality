function P = create_prob_struct( alpha, beta, n )

% "0" and "1" duration tail distribution functions and probability mass functions
P_R = [1,beta*(1-alpha).^(0:n-2)];
P_L = [1,alpha*(1-beta).^(0:n-2)];
p_R = [1-beta,beta*alpha*(1-alpha).^(0:n-2)];
p_L = [1-alpha,beta*alpha*(1-beta).^(0:n-2)];

L  = create_mtx( n, P_R, p_R ); % Pr(m-1 losses occur in the next n-1 following a loss)
L0 = create_mtx( n, p_R, p_R ); % Pr(m-1 losses occur in the next n-1 between 2 losses)
L1 = L-L0; % Pr(m-1 losses occur in the next n-1 between a loss and a reception)

R  = create_mtx( n, P_L, p_L ); % Pr(m-1 receptions occur in the next n-1 following a reception)
R1 = create_mtx( n, p_L, p_L ); % Pr(m-1 receptions occur in the next n-1 between 2 receptions)
R0 = R-R1; % Pr(m-1 receptions occur in the next n-1 between a reception and a loss)

T = [1-beta,   alpha;
       beta, 1-alpha];
   
P = struct('L0',L0,'L1',L1,'R0',R0,'R1',R1,'T',T);

return