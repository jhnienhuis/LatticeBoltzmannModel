function g = Stream(p,g)

if p.A ~= 1 % Lax-Wendroff scheme according to Guo, Shi and Zheng 2011
    
    alpha0 = 1-p.A^2;
    alphap = p.A*(p.A-1)/2;
    alpham = p.A*(p.A+1)/2;
    
    for i=1:p.Q
        g.fIn(i,:,:) = alpham*FastCircShift(g.fOut(i,:,:), [0,g.cx(i),g.cy(i)]) ...
            + alphap*FastCircShift(g.fOut(i,:,:), [0,-g.cx(i),-g.cy(i)]) ...
            + alpha0*g.fOut(i,:,:);
    end    
    
else % Standard explicit finite difference
    fIn = g.fIn;
    for i=1:p.Q
        fIn(i,:,:) = FastCircShift(g.fOut(i,:,:), [0,g.cx(i),g.cy(i)]);
    end
    g.fIn = fIn;
    
end



