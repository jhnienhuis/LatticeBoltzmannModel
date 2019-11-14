function g = Collide(p,g)



% COLLISION STEP 

rhomthreehalvesgux2pguy2 = 1 - (3/2)*(g.ux.^2+g.uy.^2);  % precalc for speed
fEq = g.fEq; % somehow keeping this separate from g allows matlab to optimize better?
for i=1:p.Q
    cu = 3*(g.cx(i)*g.ux + g.cy(i)*g.uy); 
    %    g.fEq(i,:,:) = g.t(i)*(g.rho .* ( 1 + cu + 1/2*(cu.*cu) - 3/2*(gux2pguy2) )); 
    fEq(i,:,:) = g.t(i)*(cu + 1/2*(cu.*cu) + rhomthreehalvesgux2pguy2); 
    %g.force(i,:,:) = (3*g.t(i)*(g.cx(i)*g.F(1)+g.cy(i)*g.F(2))) * g.rho;  % scalar * array
    %  g.force(i,:,:) = 3*g.t(i)*(g.cx(i)*g.F(1)+g.cy(i)*g.F(2));  % scalar * array
end
g.fEq = fEq;
g.force = 3*g.t(:).*(g.F(1)*g.cx(:) + g.F(2).*g.cy(:));
g.fNeq = g.fIn - g.fEq; % used by GetLESOmega

if p.doLES % do large-eddy simulation
    g.omega = GetLESOmega(p,g);
end
    
if ~isfield(p,'MRT') || p.MRT == 0 % no MRT
    
    %g.fOut  = g.fIn - repmat(g.omega,[p.Q,1,1]) .* (g.fIn-g.fEq) + g.force;
    %g.fOut  = g.fIn - bsxfun( @times, g.omega, g.fIn-g.fEq ) + g.force;
    g.fOut  = g.fIn - bsxfun( @times, g.omega, g.fNeq );
    g.fOut = bsxfun( @plus, g.fOut, g.force );

else % do MRT
       
    g.rIn = g.M*reshape( g.fIn, 9, p.ly*p.lx );
    % r = (rho, e, epsilon, jx, qx, jy, qy, pxx, pxy)
    % see Lallemand and Luo, JCP 2003
    jx2 = g.rIn(4,:).^2;
    jy2 = g.rIn(6,:).^2;
    %jx2jy2rho = (jx2 + jy2)./g.rIn(1,:);
    jx2jy2rho = (jx2 + jy2); % incompressible
    g.rEq(2,:) = -2*g.rIn(1,:) + 3*jx2jy2rho;
    g.rEq(3,:) = g.rIn(1,:)    - 3*jx2jy2rho;
    g.rEq(5,:) = -g.rIn(4,:);
    g.rEq(7,:) = -g.rIn(6,:);
    % g.rEq(8,:) = (jx2-jy2)./g.rIn(1,:);
    % g.rEq(9,:) = g.rIn(4,:).*g.rIn(6,:)./g.rIn(1,:);
    g.rEq(8,:) = (jx2-jy2);
    g.rEq(9,:) = g.rIn(4,:).*g.rIn(6,:);

    g.s89 = reshape(g.omega, 1, p.ly*p.lx);
    

    g.rNeq = g.rIn - g.rEq;
    
    g.rOut(1,:) = g.rIn(1,:);
    g.rOut(2,:) = g.rIn(2,:) - g.s2.*g.rNeq(2,:);
    g.rOut(3,:) = g.rIn(3,:) - g.s3.*g.rNeq(3,:);
    g.rOut(4,:) = g.rIn(4,:);
    g.rOut(5,:) = g.rIn(5,:) - g.s57.*g.rNeq(5,:);
    g.rOut(6,:) = g.rIn(6,:);
    g.rOut(7,:) = g.rIn(7,:) - g.s57.*g.rNeq(7,:);
    g.rOut(8,:) = g.rIn(8,:) - g.s89.*g.rNeq(8,:);
    g.rOut(9,:) = g.rIn(9,:) - g.s89.*g.rNeq(9,:);

    g.fOut = reshape( g.M\g.rOut, 9, p.lx, p.ly );
    g.fOut = bsxfun( @plus, g.fOut, g.force );
 
end