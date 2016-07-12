function A = create_mtx( n, B, b )
% Make sure B and b have correct length (row vectors of size(A,1))
% B: row vector st B(i)=Pr(0^(i-1) | 1) or Pr(1^(i-1) | 0)
% b: row vector st b(i)=Pr(0^(i-1)1 | 1) or Pr(1^(i-1)0 | 0)
A = zeros(n);
A(1,:) = B;
for i = 2:n
    for j = i:n
        A(i,j) = fliplr(b(1:(j-i+1))) * A(i-1,i-1:j-1)';
    end
end

return