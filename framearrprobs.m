function Phi = framearrprobs( L0, L1, R0, R1, T, k, m )
% 
N = length(k);
Phi = zeros(2,N);

for i = 1:N
    
    % find the reference frame, ref<0 iff i==1
    z = mod(i,4);
    if z == 1
        ref = i-4;
    elseif z == 2 || z == 0
        ref = i-1;
    else
        ref = i-2;
    end
    
    A = [ sum( L0(1 : m(i)-1, k(i)+m(i)-1) ), sum( R0(k(i) : k(i)+m(i)-1, k(i)+m(i)-1) );
          sum( L1(1 : m(i), k(i)+m(i)-1) ), sum( R1(k(i)-1 : k(i)+m(i)-1, k(i)+m(i)-1) ) ];
      
    if ref < 0;
        Phi(:,1) = A * [T(1,2);T(2,1)]/sum(T(1,2)+T(2,1));
    else
        inb = ref:i; inb(1) = []; inb(end) = [];
        Phi(:,i) = A * T^(sum(k(inb)+m(inb))+1) * Phi(:,ref)/sum(Phi(:,ref));
    end
    
end
Phi = (sum(Phi))';
return