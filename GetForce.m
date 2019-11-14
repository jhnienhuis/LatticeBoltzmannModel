function g = GetForce(p,g)

g.F = p.F(g.n+1,:)/p.f0; % n+1 because the first entry is for the initial state
