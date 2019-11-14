function g = MacroBC(p,g)
%According to Zhang & Kwok Phys Rev E 2006

% left boundary
in = 1;

% right boundary
out = p.lx;

% correct with average density across left and right boundary
rho_in = mean(g.rho(in,:));
rho_out= mean(g.rho(out,:));


fact_in = (1+(3*p.beta))/rho_out;
fact_out = (1-(3*p.lx*p.x0*p.beta))/rho_in;

iy = 2:(p.ly);

% horizontal velocities across boundary
g.fIn(4,in,iy) = g.fOut(2,out,iy).*fact_in;
g.fIn(2,out,iy) = g.fOut(4,in,iy).*fact_out;

% diagonal velocities across boundary
g.fIn(7,in,iy(1:end-1)) = g.fOut(9,out,iy(2:end)).*fact_in;
g.fIn(9,out,iy(2:end)) = g.fOut(7,in,iy(1:end-1)).*fact_out;

g.fIn(8,in,iy(2:end)) = g.fOut(6,out,iy(1:end-1)).*fact_in;
g.fIn(6,out,iy(1:end-1)) = g.fOut(8,in,iy(2:end)).*fact_out;

%simple bounce back?
%g.fOut(:,g.bbRegion) = g.fIn(g.opp,g.bbRegion); 

%g.fIn(8,in,2) = g.fOut(8,in,2);
%g.fIn(9,out,2)= g.fOut(9,out,2);

g.fIn(8,in,2)= g.fOut(9,out,2).*fact_in;
g.fIn(9,out,2)= g.fOut(8,in,2).*fact_out;


end