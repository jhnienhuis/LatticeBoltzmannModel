
function [varargout] = TraceParticle(out)

%x = zeros(1,numel(out.time));
%y = zeros(1,numel(out.time));

%for i=1:numel(out.time)
[x y] = TraceParticle2(out.velocity.ux,out.velocity.uy,out.p.lx/2,out.p.ly-10,out.p.t0*out.p.saveint);

plot(x,y);

%saveas(gcf,[out.p.savename 'traceparticle.fig'])

%close gcf

% 
% if nargout,
%     varargout = {x,y};
%     
% else,
%    plot(x,y)
%    xlabel('x (m)');
%    ylabel('y (m)');
%     
% end
    
%x&y are given in meters!!
%end


end





function [x y] = TraceParticle2(ux,uy,xi,yi,dt)

% Records the position of a tracer particle starting at position xi,yi in a
% velocity field given by components ux and uy. 
% ux and uy have dimensions [Nx Ny N], where N is the number of time slices.
% xi and yi are scalars, given in the units of Nx and Ny.
% dt is the time increment.

if length(xi) > 1
  for i = 1:length(xi)
    fprintf(1,'.');
    [x(:,i) y(:,i)] = TraceParticle2(ux,uy,xi(i),yi(i),dt);
  end
  return
end

[Nx Ny N] = size(ux);

[Y X] = meshgrid(0:Ny+1,0:Nx+1);

ux(isnan(ux)) = 0;
uy(isnan(uy)) = 0;

% xp = xi; % xp and yp are the coordinates after applying periodic boundary conditions. They will always lie on the intervals [1 Nx] and [1 Ny]
% yp = yi;

x = zeros(N,1);
y = zeros(N,1);

for n = 1:N
    x(n) = xi;
    y(n) = yi;

    [uxi uyi] = Getu(X,Y,ux(:,:,n),uy(:,:,n),xi,yi);
    
    xi = xi + dt*uxi;
    yi = yi + dt*uyi;
        
end
end
% Subfunctions

function [uxi uyi] = Getu(X,Y,ux,uy,xi,yi)

[Nx Ny N] = size(ux);
[xp yp] = PeriodicBC(xi,yi,Nx,Ny);
ux = padarray(ux,[1 1],'circular');
uy = padarray(uy,[1 1],'circular');

uxi = interp2(Y,X,ux,yp,xp);
uyi = interp2(Y,X,uy,yp,xp);

end
function [xp yp] = PeriodicBC(xi,yi,Nx,Ny)

% if xi > Nx
%     xp = xi - Nx;
% elseif xi < 1
%     xp = xi + Nx;
% else
%     xp = xi;
% end
% 
% if yi > Ny
%     yp = yi - Ny;
% elseif yi < 1
%     yp = yi + Ny;
% else
%     yp = yi;
% end

xp = mod(xi-Nx,Nx);
yp = mod(yi-Ny,Ny);
end