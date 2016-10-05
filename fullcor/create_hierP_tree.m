function t = create_hierP_tree( L, N )

gop = 2^(L-1);
G = N / gop;
t = tree(1);
for b = 2:G
    t = t.addnode(b-1,1+(b-1)*gop);
end

while nnodes(t) ~= N
    gop = gop / 2;
    for node = t.nodeorderiterator
        index = t.get( node );
        t = t.addnode( node, index + gop );
    end
end

return