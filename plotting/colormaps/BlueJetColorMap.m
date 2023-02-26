function map = BlueJetColorMap(n)

if nargin <1, n  = 100; end

n  =  round([15 40 35 30 25]./100*n);

map = makeColorMap([0  0 1],n(2),...
                   [0  1 .35],n(3),[1  1 0 ],n(4),...
                   [1 0 0],n(5), [.5 0 0 ]);
